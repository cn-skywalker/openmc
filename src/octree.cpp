#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include <openmc/bounding_box.h>
#include <openmc/cell.h>
#include <openmc/constants.h>
#include <openmc/octree.h>

namespace openmc {
bool OctreeNode::insert(const int32_t& sphere_token)
{
  // 检查球体是否在节点边界内
  if (!boundary_.intersects(sphere_token)) {
    return false;
  }

  // 如果节点未满，直接插入
  if (!divided_ && spheres_indexs_.size() < capacity_) {
    spheres_indexs_.push_back(sphere_token);
    return true;
  }

  // 如果节点已满且未分割，先分割
  if (!divided_) {
    // 如果尺寸已经小于最小分割尺寸，不再分割，直接插入当前节点
    if (!shouldSubdivide()) {
      spheres_indexs_.push_back(sphere_token);
      return true;
    }
    subdivide();

    // 重要：将当前节点的球体重新分配到子节点
    auto old_spheres = std::move(spheres_indexs_);
    spheres_indexs_.clear();

    for (const auto& token : old_spheres) {
      // 尝试插入到子节点
      bool inserted = false;
      for (auto& child : children_) {
        if (child->insert(token)) {
          inserted = true;
        }
      }
      // 如果球体无法插入任何子节点，则进行报错
      if (!inserted) {
        fatal_error("Error in OctreeNode::insert: could not reinsert existing "
                    "sphere {} into child nodes.",
          model::surfaces[abs(token) - 1]->id_);
      }
    }
  }

  // 尝试将新球体插入到子节点
  bool inserted_to_child = false;
  for (auto& child : children_) {
    if (child->insert(sphere_token)) {
      inserted_to_child = true;
    }
  }

  // 如果球体无法插入任何子节点，留在当前节点
  if (!inserted_to_child) {
    spheres_indexs_.push_back(sphere_token);
  }

  return true;
}

int32_t OctreeNode::queryPoint(const Position& point) const
{
  // 如果点不在节点边界内，返回-1
  if (!boundary_.contains(point)) {
    return -1;
  }

  // 检查当前节点的球体
  for (const auto& sphere_index : spheres_indexs_) {
    vector<double> triso_center =
      model::surfaces[abs(sphere_index) - 1]->get_center();
    double triso_radius = model::surfaces[abs(sphere_index) - 1]->get_radius();
    if (pow(point.x - triso_center[0], 2) + pow(point.y - triso_center[1], 2) +
          pow(point.z - triso_center[2], 2) <
        pow(triso_radius, 2)) {
      return sphere_index;
    }
  }

  // 检查子节点
  if (divided_) {
    for (const auto& child : children_) {
      int result = child->queryPoint(point);
      if (result != -1) {
        return result;
      }
    }
  }

  return -1;
}

std::pair<int32_t, double> OctreeNode::queryRay(
  const Position& origin, const Position& direction, int32_t on_surface) const
{
  Position current_position = origin; // 步骤（1）：令current_position=origin
  double minT = std::numeric_limits<double>::max();
  int32_t result_sphere = -1;

  // 设置最大迭代次数防止无限循环
  const int max_iterations = 1000;
  int iteration_count = 0;

  while (iteration_count++ < max_iterations) {
    // 步骤（2）：根据current_position定位当前叶子节点
    const OctreeNode* leaf_node = findLeafNode(current_position);
    if (!leaf_node) {
      // 如果不存在叶子节点，向前移动一个很小的距离再次判断
      current_position = current_position + direction * FP_PRECISION;
      leaf_node = findLeafNode(current_position);
    }

    if (!leaf_node) {
      // 如果不存在叶子节点，返回无穷大
      return {-1, std::numeric_limits<double>::max()};
    }

    // 检测该节点中所有的颗粒是否相交
    bool found_intersection = false;
    for (const auto& sphere_token : leaf_node->spheres_indexs_) {
      bool coincident = (std::abs(sphere_token) == std::abs(on_surface));
      double t = model::surfaces[abs(sphere_token) - 1]->distance(
        origin, direction, coincident);

      if (t > 0 && t < minT) {
        minT = t;
        result_sphere = sphere_token;
        found_intersection = true;
      }
    }

    // 步骤（3）：如果有颗粒相交则返回
    if (found_intersection) {
      return {result_sphere, minT};
    }

    // 步骤（4）：计算沿着射线方向从current_position到当前叶子节点出口的距离
    double exit_distance =
      leaf_node->getExitDistance(current_position, direction);

    // 如果没有出口，向前移动一小段距离，再次定位叶子节点
    if (exit_distance >= std::numeric_limits<double>::max()) {
      current_position = current_position + direction * FP_PRECISION;
      leaf_node = findLeafNode(current_position);
      if(!leaf_node){
        return {-1, std::numeric_limits<double>::max()};
      }
      exit_distance = leaf_node->getExitDistance(current_position, direction);
    }

    if (exit_distance >= std::numeric_limits<double>::max()) {
      // 没有出口，射线离开八叉树边界
      return {-1, std::numeric_limits<double>::max()};
    }

    // 将current_position移动到出口位置
    current_position =
      current_position + direction * (exit_distance + FP_COINCIDENT);

    // 步骤（5）：回到步骤2（通过循环实现）
  }

  // 达到最大迭代次数，返回无交点
  return {-1, std::numeric_limits<double>::max()};
}

std::pair<int32_t, double> OctreeNode::queryRayold(
  const Position& origin, const Position& direction, int32_t on_surface) const
{
  // 如果射线与节点边界不相交，返回无效结果
  if (!boundary_.rayIntersect(origin, direction)) {
    return {-1, std::numeric_limits<double>::max()};
  }

  std::pair<int32_t, double> result = {-1, std::numeric_limits<double>::max()};
  double minT = std::numeric_limits<double>::max();

  // 检查当前节点的球体
  for (const auto& sphere_token : spheres_indexs_) {
    bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
    double t = model::surfaces[abs(sphere_token) - 1]->distance(
      origin, direction, coincident);
    if (t > 0 && t < minT) {
      minT = t;
      result.first = sphere_token;
      result.second = t;
    }
  }

  // 检查子节点
  if (divided_) {
    for (const auto& child : children_) {
      auto childResult = child->queryRay(origin, direction, on_surface);
      if (childResult.first != -1 && childResult.second < minT) {
        minT = childResult.second;
        result = childResult;
      }
    }
  }

  return result;
}

void OctreeNode::subdivide()
{
  Position center = boundary_.getCenter();
  Position min = boundary_.min();
  Position max = boundary_.max();

  children_.reserve(8);
  int child_depth = depth_ + 1;

  // 统一使用min/max/center来定义边界，确保无重叠无遗漏
  // 为每个子节点计算Morton编码
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, min.z),
                                   Position(center.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(0), child_depth));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(center.x, min.y, min.z),
                                   Position(max.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(1), child_depth));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, center.y, min.z),
                                   Position(center.x, max.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(2), child_depth));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, min.z), Position(max.x, max.y, center.z)),
    min_size_, capacity_, computeChildMortonCode(3), child_depth));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, center.z),
                                   Position(center.x, center.y, max.z)),
      min_size_, capacity_, computeChildMortonCode(4), child_depth));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, min.y, center.z), Position(max.x, center.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(5), child_depth));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(min.x, center.y, center.z), Position(center.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(6), child_depth));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, center.z), Position(max.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(7), child_depth));

  divided_ = true;
}

void OctreeNode::printTree(int depth, bool showAll) const
{
  if (!showAll && spheres_indexs_.empty() && !divided_) {
    return; // 如果节点为空且未分割，且不要求显示所有节点，则跳过
  }
  // 创建缩进字符串
  std::string indent(depth * 2, ' ');
  // 打印当前节点信息
  std::cout << indent << "└─ Node [depth=" << depth
            << ", spheres=" << spheres_indexs_.size()
            << ", divided=" << (divided_ ? "true" : "false") << "]\n";
  if (showAll) {

    // 打印边界框信息
    Position center = boundary_.getCenter();
    Position size = boundary_.getSize();
    std::cout << indent << "   Bounds: min(" << boundary_.xmin << ", "
              << boundary_.ymin << ", " << boundary_.zmin << ")"
              << " max(" << boundary_.xmax << ", " << boundary_.ymax << ", "
              << boundary_.zmax << ")\n";
    std::cout << indent << "   Center: (" << center.x << ", " << center.y
              << ", " << center.z << ")"
              << " Size: (" << size.x << ", " << size.y << ", " << size.z
              << ")\n";
  }

  // 打印当前节点中的球体
  if (!spheres_indexs_.empty()) {
    std::cout << indent << "   Spheres: ";
    for (size_t i = 0; i < spheres_indexs_.size(); ++i) {
      std::cout << spheres_indexs_[i];
      if (i < spheres_indexs_.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "\n";
  }

  // 递归打印子节点
  if (divided_) {
    for (size_t i = 0; i < children_.size(); ++i) {
      std::cout << indent << "  Child " << i + 1 << ":\n";
      children_[i]->printTree(depth + 1);
    }
  }
}

const OctreeNode* OctreeNode::findLeafNode(const Position& point) const
{
  // 如果点不在节点边界内，返回 nullptr
  if (!boundary_.contains(point)) {
    return nullptr;
  }

  // 如果当前节点是叶子节点（未分割），返回当前节点
  if (!divided_) {
    return this;
  }

  // 如果节点已分割，递归查询子节点
  for (const auto& child : children_) {
    const OctreeNode* result = child->findLeafNode(point);
    if (result != nullptr) {
      return result;
    }
  }

  // 理论上不会执行到这里，因为点应该在某个子节点中
  return nullptr;
}

double OctreeNode::getExitDistance(
  const Position& origin, const Position& direction) const
{
  auto [t_enter, t_exit] =
    boundary_.rayIntersectionDistances(origin, direction);

  if (t_exit < std::numeric_limits<double>::max()) {
    return t_exit;
  }

  return std::numeric_limits<double>::max();
}

uint64_t OctreeNode::computeChildMortonCode(int child_index) const
{
  // 子节点的Morton编码 = 父节点编码左移3位 + 子节点索引
  return (morton_code_ << 3) | (child_index & 0x7);
}

bool OctreeNode::shouldSubdivide() const
{
  // 检查边界尺寸是否大于最小分割尺寸
  Position size = boundary_.getSize();
  return (size.x > min_size_ && size.y > min_size_ && size.z > min_size_);
}

} // namespace openmc