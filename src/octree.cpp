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
  // 如果射线与节点边界不相交，返回无效结果-1
  if (!rayAABBIntersect(origin, direction, boundary_)) {
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

  // 统一使用min/max/center来定义边界，确保无重叠无遗漏
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, min.z),
                                   Position(center.x, center.y, center.z)),
      capacity_));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(center.x, min.y, min.z),
                                   Position(max.x, center.y, center.z)),
      capacity_));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, center.y, min.z),
                                   Position(center.x, max.y, center.z)),
      capacity_));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, min.z), Position(max.x, max.y, center.z)),
    capacity_));

  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, center.z),
                                   Position(center.x, center.y, max.z)),
      capacity_));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, min.y, center.z), Position(max.x, center.y, max.z)),
    capacity_));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(min.x, center.y, center.z), Position(center.x, max.y, max.z)),
    capacity_));

  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, center.z), Position(max.x, max.y, max.z)),
    capacity_));

  divided_ = true;
}
bool OctreeNode::rayAABBIntersect(const Position& origin,
  const Position& direction, const BoundingBox& aabb) const
{
  double tmin = 0.0; // 射线 t >= 0
  double tmax = std::numeric_limits<double>::max();

  for (int i = 0; i < 3; ++i) {
    double dir = direction[i];
    double minVal = aabb.min()[i];
    double maxVal = aabb.max()[i];

    if (std::abs(dir) < FP_PRECISION) {
      // 射线平行于该轴
      if (origin[i] < minVal || origin[i] > maxVal) {
        return false; // 不在范围内，不相交
      }
      // 否则不限制 t 范围，跳过
    } else {
      double invD = 1.0 / dir;
      double t0 = (minVal - origin[i]) * invD;
      double t1 = (maxVal - origin[i]) * invD;

      if (invD < 0.0) {
        std::swap(t0, t1);
      }

      tmin = std::max(t0, tmin);
      tmax = std::min(t1, tmax);

      if (tmax <= tmin) {
        return false;
      }
    }
  }

  return true;
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

} // namespace openmc