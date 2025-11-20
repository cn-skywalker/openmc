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
namespace model {
std::unordered_map<uint64_t, OctreeNode*> leaf_nodes_map;
} // namespace model

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
    fatal_error("Error in OctreeNode::insert: could not reinsert existing "
                "sphere {} into child nodes.",
      model::surfaces[abs(sphere_token) - 1]->id_);
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

std::pair<int32_t, double> OctreeNode::queryRay_morton_code(
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // 步骤（1）：令current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  uint64_t old_morton_code = -1;
  bool if_stuck = false;
  int32_t stuck_count = 0;
  BoxFace exit_face = BoxFace::NONE;
  double total_distance = 0.0;
  const OctreeNode* leaf_node = nullptr;

  while (true) {
    bool found = false; // 标记是否被找到了
    // 先用morton码进行查找
    if (old_morton_code != -1 && leaf_node != nullptr) {
      auto [next_morton_code, next_depth] =
        leaf_node->deriveNextNodeMorton(exit_face);

      // 对深度进行循环回退，直到找到对应的叶子节点或无法回退为止
      for (int depth = next_depth; depth > 0; depth--) {
        uint64_t adjusted_code = next_morton_code >> (3 * (next_depth - depth));
        auto it = model::leaf_nodes_map.find(adjusted_code);
        if (it != model::leaf_nodes_map.end()) {
          leaf_node = it->second;
          if (leaf_node->divided_ == false) {
            found = true;
          }
          break;
        }
      }
    }
    if (!found) {
      // 如果没找到，使用位置查找
      // 查找当前位置所在的叶子节点
      if (leaf_node == nullptr) {
        leaf_node = findLeafNode(current_position);
      } else {
        leaf_node = leaf_node->findLeafNode(current_position);
      }
      if (leaf_node == nullptr ||
          leaf_node->getMortonCode() == old_morton_code) {
        // 如果不在任何叶子节点内，向前移动一小段距离
        current_position += direction * FP_COINCIDENT;
        // 再次全局查找，如果还是不在则退出
        leaf_node = findLeafNode(current_position);
        if (leaf_node == nullptr) {
          break;
        }
        if_stuck = leaf_node->getMortonCode() == old_morton_code;
        if (if_stuck) {
          stuck_count++;
          if (stuck_count > 10) {
            fatal_error("Fatal error in OctreeNode::queryRay: stuck in the "
                        "same leaf node "
                        "too many times.");
          }
        }
      }
    }
    // 如果没有堵在同一个叶子节点，查找当前位置的叶子节点内的球体
    if (!if_stuck) {
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // 记录球体token
        }
      }
    }
    // 计算射线与当前节点边界的出口距离
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // 如果出口距离无限大，说明射线不会再进入其他节点，退出循环
    exit_face = exit_face1;
    if (exit_distance == INFTY) {
      break;
    }
    // 步骤（5）：更新current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    total_distance += exit_distance + FP_COINCIDENT;
    // 如果超过最大距离，退出循环
    if (total_distance > max_distance) {
      break;
    }
    // 更新old_morton_code
    old_morton_code = leaf_node->getMortonCode();
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRay_Neighbor_search(
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // 步骤（1）：令current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  double total_distance = 0.0;
  BoxFace exit_face = BoxFace::NONE;

  const OctreeNode* leaf_node = nullptr; // 当前位置所在的叶子节点

  while (true) {
    bool found = false;
    bool if_stuck = false;
    // 如果出口不为空，且存在旧节点
    if (exit_face != BoxFace::NONE && old_leaf_node != nullptr) {
      auto& Neighbor_nodes = old_leaf_node->getNeighbors(exit_face);
      // 如果邻居列表不为空
      if (!Neighbor_nodes.empty()) {
        // 判断是否有邻居节点包含current_position
        for (auto neighbor : Neighbor_nodes) {
          if (neighbor->boundary_.contains(current_position)) {
            leaf_node = neighbor;
            found = true;
            break;
          }
        }
      }
    }

    if (!found) {
      // 查找当前位置所在的叶子节点
      leaf_node = findLeafNode(current_position);
      if (leaf_node == nullptr || leaf_node == old_leaf_node) {
        // 如果不在任何叶子节点内，向前移动一小段距离
        current_position += direction * FP_COINCIDENT;
        // 再次查找，如果还是不在则退出
        leaf_node = findLeafNode(current_position);
        if (leaf_node == nullptr) {
          break;
        }
        if_stuck = leaf_node == old_leaf_node;
        if (if_stuck) {
          warning(
            "Warning in OctreeNode::queryRay: stuck in the same leaf node.");
        }
      } else {
        // 如果没有堵塞，添加新节点到旧节点的邻居列表
        if (old_leaf_node != nullptr) {
          old_leaf_node->addNeighbor(
            exit_face, const_cast<OctreeNode*>(leaf_node));
        }
      }
    }

    // 如果没有堵在同一个叶子节点，查找当前位置的叶子节点内的球体
    if (!if_stuck) {
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // 记录球体token
        }
      }
    }
    // 计算射线与当前节点边界的出口距离
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    exit_face = exit_face1;
    // 如果出口距离无限大，说明射线不会再进入其他节点，退出循环

    if (exit_distance == INFTY) {
      break;
    }
    // 步骤（5）：更新current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    total_distance += exit_distance + FP_COINCIDENT;
    // 如果超过最大距离，退出循环
    if (total_distance > max_distance) {
      break;
    }
    old_leaf_node = const_cast<OctreeNode*>(leaf_node);
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRay_leaf_find(
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // 步骤（1）：令current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  bool if_stuck = false;
  int32_t stuck_count = 0;

  while (true) {
    // 查找当前位置所在的叶子节点
    const OctreeNode* leaf_node = findLeafNode(current_position);
    if (leaf_node == nullptr || leaf_node == old_leaf_node) {
      // 如果不在任何叶子节点内，向前移动一小段距离
      current_position += direction * FP_COINCIDENT;
      // 再次查找，如果还是不在则退出
      leaf_node = findLeafNode(current_position);
      if (leaf_node == nullptr) {
        break;
      }
      if_stuck = leaf_node == old_leaf_node;
      if (if_stuck) {
        stuck_count++;
        if (stuck_count > 10) {
          fatal_error(
            "Fatal error in OctreeNode::queryRay: stuck in the same leaf node "
            "too many times.");
        }
      }
    }
    // 如果没有堵在同一个叶子节点，查找当前位置的叶子节点内的球体
    if (!if_stuck) {
      stuck_count = 0;
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // 记录球体token
        }
      }
    }
    // 计算射线与当前节点边界的出口距离
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // 如果出口距离无限大，说明射线不会再进入其他节点，退出循环
    if (exit_distance == INFTY) {
      break;
    }
    // 步骤（5）：更新current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    old_leaf_node = const_cast<OctreeNode*>(leaf_node);
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRay_leaf_find_old(
  // 没有定位加速的
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // 步骤（1）：令current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  bool if_stuck = false;
  int32_t stuck_count = 0;

  while (true) {
    // 查找当前位置所在的叶子节点
    const OctreeNode* leaf_node = findLeafNode_old(current_position);
    if (leaf_node == nullptr || leaf_node == old_leaf_node) {
      // 如果不在任何叶子节点内，向前移动一小段距离
      current_position += direction * FP_COINCIDENT;
      // 再次查找，如果还是不在则退出
      leaf_node = findLeafNode_old(current_position);
      if (leaf_node == nullptr) {
        break;
      }
      if_stuck = leaf_node == old_leaf_node;
      if (if_stuck) {
        stuck_count++;
        if (stuck_count > 10) {
          fatal_error(
            "Fatal error in OctreeNode::queryRay: stuck in the same leaf node "
            "too many times.");
        }
      }
    }
    // 如果没有堵在同一个叶子节点，查找当前位置的叶子节点内的球体
    if (!if_stuck) {
      stuck_count = 0;
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // 记录球体token
        }
      }
    }
    // 计算射线与当前节点边界的出口距离
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // 如果出口距离无限大，说明射线不会再进入其他节点，退出循环
    if (exit_distance == INFTY) {
      break;
    }
    // 步骤（5）：更新current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    old_leaf_node = const_cast<OctreeNode*>(leaf_node);
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRayold(
  const Position& origin, const Position& direction, int32_t on_surface) const
{
  // 如果射线与节点边界不相交，返回无效结果
  if (!boundary_.rayIntersect(origin, direction)) {
    return {std::numeric_limits<int32_t>::max(), std::numeric_limits<double>::max()};
  }

  std::pair<int32_t, double> result = {std::numeric_limits<int32_t>::max(), std::numeric_limits<double>::max()};
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
      auto childResult = child->queryRayold(origin, direction, on_surface);
      if (childResult.first != std::numeric_limits<int32_t>::max() && childResult.second < minT) {
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
  // 为每个子节点计算Morton编码，右手系建模
  // 正向为1，负向为0，编码顺序为zyx

  // 索引 0: 左-下-前 (000)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, min.z),
                                   Position(center.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(0), child_depth));

  // 索引 1: 右-下-前 (001)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(center.x, min.y, min.z),
                                   Position(max.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(1), child_depth));

  // 索引 2: 左-上-前 (010)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, center.y, min.z),
                                   Position(center.x, max.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(2), child_depth));

  // 索引 3: 右-上-前 (011)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, min.z), Position(max.x, max.y, center.z)),
    min_size_, capacity_, computeChildMortonCode(3), child_depth));

  // 索引 4: 左-下-后 (100)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, center.z),
                                   Position(center.x, center.y, max.z)),
      min_size_, capacity_, computeChildMortonCode(4), child_depth));

  // 索引 5: 右-下-后 (101)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, min.y, center.z), Position(max.x, center.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(5), child_depth));

  // 索引 6: 左-上-后 (110)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(min.x, center.y, center.z), Position(center.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(6), child_depth));

  // 索引 7: 右-上-后 (111)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, center.z), Position(max.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(7), child_depth));

  divided_ = true;
}

void OctreeNode::printBasicInfo() const
{
  // 打印当前树结构的最深深度和叶子节点数量
  int max_depth = 0;
  std::function<void(const OctreeNode*, int)> traverse;
  traverse = [&](const OctreeNode* node, int depth) {
    if (!node->divided_) {

      if (depth > max_depth) {
        max_depth = depth;
      }
    } else {
      for (const auto& child : node->children_) {
        traverse(child.get(), depth + 1);
      }
    }
  };
  traverse(this, 0);

  // 输出最深深度
  write_message(fmt::format(
    "Octree basic info: max depth={}, max capacity={}", max_depth, capacity_));
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

std::vector<const OctreeNode*> OctreeNode::findLeafNodesContainingSphere(
  int32_t sphere_token) const
{
  std::vector<const OctreeNode*> result;

  // 检查当前节点是否包含该球体
  bool contains_sphere = false;
  for (const auto& token : spheres_indexs_) {
    if (token == sphere_token) {
      contains_sphere = true;
      break;
    }
  }

  // 如果当前节点是叶子节点且包含该球体，添加到结果
  if (!divided_ && contains_sphere) {
    result.push_back(this);
  }

  // 如果节点已分割，递归检查子节点
  if (divided_) {
    for (const auto& child : children_) {
      auto child_results = child->findLeafNodesContainingSphere(sphere_token);
      result.insert(result.end(), child_results.begin(), child_results.end());
    }
  }
  std::cout << "Sphere " << sphere_token << " is contained in ";

  std::cout << "Sphere " << sphere_token << " is contained in " << result.size()
            << " leaf nodes:\n";
  for (const auto* node : result) {
    Position center = node->boundary_.getCenter();
    Position size = node->boundary_.getSize();
    std::cout << "  Leaf node at (" << center.x << ", " << center.y << ", "
              << center.z << ") size (" << size.x << ", " << size.y << ", "
              << size.z << ") depth=" << node->getDepth() << "\n";
  }

  return result;
}

const OctreeNode* OctreeNode::findLeafNode(const Position& point) const
{
  // 如果点不在节点边界内，返回 nullptr
  if (!boundary_.contains(point)) {
    return nullptr;
  }

  const OctreeNode* current = this;

  // 迭代向下查找，直到叶子节点
  while (current->divided_) {
    Position center = current->boundary_.getCenter();

    // 根据点相对于中心的位置确定子节点索引
    int child_index = 0;
    if (point.x >= center.x)
      child_index |= 1; // x 轴：0=左, 1=右
    if (point.y >= center.y)
      child_index |= 2; // y 轴：0=下, 1=上
    if (point.z >= center.z)
      child_index |= 4; // z 轴：0=前, 1=后

    // 直接访问对应的子节点
    current = current->children_[child_index].get();

    // 安全检查：如果子节点为空，返回当前节点
    if (current == nullptr) {
      return const_cast<OctreeNode*>(this);
    }
  }

  return current;
}

const OctreeNode* OctreeNode::findLeafNode_old(const Position& point) const
{
  // 如果点不在节点边界内，返回 nullptr
  if (!boundary_.contains(point)) {
    return nullptr;
  }

  const OctreeNode* current = this;

  // 迭代向下查找，直到叶子节点
  while (current->divided_) {
    // 暴力判断：逐个检查每个子节点
    bool found = false;

    // 检查所有8个子节点，找到包含该点的那个
    for (int i = 0; i < 8; i++) {
      if (current->children_[i] != nullptr &&
          current->children_[i]->boundary_.contains(point)) {
        current = current->children_[i].get();
        found = true;
        break;
      }
    }

    // 如果没有找到包含该点的子节点，返回当前节点
    if (!found) {
      return current;
    }
  }

  return current;
}

std::pair<BoxFace, double> OctreeNode::getExitDistance(
  const Position& origin, const Position& direction) const
{
  return boundary_.rayIntersectionDistances(origin, direction);
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

std::pair<uint64_t, int> OctreeNode::deriveNextNodeMorton(
  BoxFace exit_face) const
{
  uint64_t next_morton_code = morton_code_;

  // 根据出口面直接调整 Morton 码
  switch (exit_face) {
  case BoxFace::MIN_X:
    // X坐标减1：需要找到最低的x位为1的位置，将其变为0，并将所有更低的x位变为1
    for (int i = 0; i < depth_; ++i) {
      uint64_t x_bit_pos = 3 * i;
      if ((next_morton_code >> x_bit_pos) & 0x1) {
        // 找到第一个为1的x位，将其变为0
        next_morton_code &= ~(1ULL << x_bit_pos);
        // 将所有更低的x位设为1
        for (int j = 0; j < i; ++j) {
          next_morton_code |= (1ULL << (3 * j));
        }
        break;
      }
    }
    break;

  case BoxFace::MAX_X:
    // X坐标加1：需要找到最低的x位为0的位置，将其变为1，并将所有更低的x位变为0
    for (int i = 0; i < depth_; ++i) {
      uint64_t x_bit_pos = 3 * i;
      if (!((next_morton_code >> x_bit_pos) & 0x1)) {
        // 找到第一个为0的x位，将其变为1
        next_morton_code |= (1ULL << x_bit_pos);
        // 将所有更低的x位设为0
        for (int j = 0; j < i; ++j) {
          next_morton_code &= ~(1ULL << (3 * j));
        }
        break;
      }
    }
    break;

  case BoxFace::MIN_Y:
    // Y坐标减1：操作y位（位置 3*i+1）
    for (int i = 0; i < depth_; ++i) {
      uint64_t y_bit_pos = 3 * i + 1;
      if ((next_morton_code >> y_bit_pos) & 0x1) {
        next_morton_code &= ~(1ULL << y_bit_pos);
        for (int j = 0; j < i; ++j) {
          next_morton_code |= (1ULL << (3 * j + 1));
        }
        break;
      }
    }
    break;

  case BoxFace::MAX_Y:
    // Y坐标加1：操作y位（位置 3*i+1）
    for (int i = 0; i < depth_; ++i) {
      uint64_t y_bit_pos = 3 * i + 1;
      if (!((next_morton_code >> y_bit_pos) & 0x1)) {
        next_morton_code |= (1ULL << y_bit_pos);
        for (int j = 0; j < i; ++j) {
          next_morton_code &= ~(1ULL << (3 * j + 1));
        }
        break;
      }
    }
    break;

  case BoxFace::MIN_Z:
    // Z坐标减1：操作z位（位置 3*i+2）
    for (int i = 0; i < depth_; ++i) {
      uint64_t z_bit_pos = 3 * i + 2;
      if ((next_morton_code >> z_bit_pos) & 0x1) {
        next_morton_code &= ~(1ULL << z_bit_pos);
        for (int j = 0; j < i; ++j) {
          next_morton_code |= (1ULL << (3 * j + 2));
        }
        break;
      }
    }
    break;

  case BoxFace::MAX_Z:
    // Z坐标加1：操作z位（位置 3*i+2）
    for (int i = 0; i < depth_; ++i) {
      uint64_t z_bit_pos = 3 * i + 2;
      if (!((next_morton_code >> z_bit_pos) & 0x1)) {
        next_morton_code |= (1ULL << z_bit_pos);
        for (int j = 0; j < i; ++j) {
          next_morton_code &= ~(1ULL << (3 * j + 2));
        }
        break;
      }
    }
    break;

  default:
    break;
  }

  return {next_morton_code, depth_};
}

void buildLeafMap(OctreeNode* node)
{
  model::leaf_nodes_map[node->getMortonCode()] = node;
  if (node->divided_) {
    for (auto& child : node->children_) {
      buildLeafMap(child.get());
    }
  }
}
} // namespace openmc