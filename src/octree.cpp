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
      // 如果球体无法插入任何子节点，留在当前节点
      if (!inserted) {
        spheres_indexs_.push_back(token);
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

uint64_t OctreeNode::computeChildMortonCode(int child_index) const
{
  // 子节点的Morton编码 = 父节点编码左移3位 + 子节点索引
  return (morton_code_ << 3) | (child_index & 0x7);
}

uint64_t OctreeNode::computeMortonCode(
  const Position& pos, const Position& min, const Position& max, int max_depth)
{
  // 将位置归一化到[0,1]范围
  double x_norm = (pos.x - min.x) / (max.x - min.x);
  double y_norm = (pos.y - min.y) / (max.y - min.y);
  double z_norm = (pos.z - min.z) / (max.z - min.z);

  // 将归一化坐标映射到整数空间
  uint32_t x_int = static_cast<uint32_t>(x_norm * ((1 << max_depth) - 1));
  uint32_t y_int = static_cast<uint32_t>(y_norm * ((1 << max_depth) - 1));
  uint32_t z_int = static_cast<uint32_t>(z_norm * ((1 << max_depth) - 1));

  // 计算Morton编码（交错位）
  uint64_t code = 0;
  for (int i = 0; i < max_depth; ++i) {
    code |= ((x_int >> i) & 1) << (3 * i);
    code |= ((y_int >> i) & 1) << (3 * i + 1);
    code |= ((z_int >> i) & 1) << (3 * i + 2);
  }

  return code;
}

void OctreeNode::decodeMortonCode(uint64_t code, int depth, Position& min,
  Position& max, const Position& root_min, const Position& root_max)
{
  double size_x = (root_max.x - root_min.x) / (1 << depth);
  double size_y = (root_max.y - root_min.y) / (1 << depth);
  double size_z = (root_max.z - root_min.z) / (1 << depth);

  uint64_t temp_code = code;
  int x_idx = 0, y_idx = 0, z_idx = 0;

  for (int i = 0; i < depth; ++i) {
    x_idx |= (temp_code & 1) << i;
    y_idx |= ((temp_code >> 1) & 1) << i;
    z_idx |= ((temp_code >> 2) & 1) << i;
    temp_code >>= 3;
  }

  min.x = root_min.x + x_idx * size_x;
  min.y = root_min.y + y_idx * size_y;
  min.z = root_min.z + z_idx * size_z;
  max.x = min.x + size_x;
  max.y = min.y + size_y;
  max.z = min.z + size_z;
}

bool OctreeNode::shouldSubdivide() const
{
  // 检查边界尺寸是否大于最小分割尺寸
  Position size = boundary_.getSize();
  return (size.x > min_size_ && size.y > min_size_ && size.z > min_size_);
}

} // namespace openmc