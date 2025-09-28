#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include <openmc/bounding_box.h>
#include <openmc/cell.h>
#include <openmc/octree.h>

namespace openmc {
bool OctreeNode::insert(const int32_t& sphere_token)
{
  // 如果球体不在节点边界内，不插入
  if (!boundary_.intersects(sphere_token)) {
    return false;
  }

  // 如果节点未满，直接插入
  if (spheres_indexs_.size() < capacity_) {
    spheres_indexs_.push_back(sphere_token);
    return true;
  }

  // 如果节点已满且未分割，先分割
  if (!divided_) {
    subdivide();
  }

  // 尝试将球体插入子节点
  for (auto& child : children_) {
    if (child->insert(sphere_token)) {
      return true;
    }
  }

  // 如果球体无法插入任何子节点，留在当前节点
  spheres_indexs_.push_back(sphere_token);
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

std::pair<int, Position> OctreeNode::queryRay(
  const Position& origin, const Position& direction, int32_t on_surface) const
{
  // 如果射线与节点边界不相交，返回-1
  if (!rayAABBIntersect(origin, direction, boundary_)) {
    return {-1, Position()};
  }

  std::pair<int, Position> result = {-1, Position()};
  double minT = std::numeric_limits<double>::max();

  // 检查当前节点的球体
  for (const auto& sphere_token : spheres_indexs_) {
    bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
    double t = model::surfaces[abs(sphere_token) - 1]->distance(
      origin, direction, coincident);
    if (t > 0 && t < minT) {
      minT = t;
      result.first = sphere_token;
      result.second = origin + direction * t;
    }
  }

  // 检查子节点
  if (divided_) {
    for (const auto& child : children_) {
      auto childResult = child->queryRay(origin, direction, on_surface);
      if (childResult.first != -1) {
        Position hitPoint = childResult.second;
        Position l = (hitPoint - origin);
        double t = sqrt(l.x * l.x + l.y * l.y + l.z * l.z);
        if (t < minT) {
          minT = t;
          result = childResult;
        }
      }
    }
  }

  return result;
}

void OctreeNode::subdivide()
{
  Position center = boundary_.getCenter();
  Position size = boundary_.getSize();
  Position halfSize = size * 0.5;

  // 创建8个子节点
  children_.reserve(8);

  // 前下左
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(boundary_.min(), center), capacity_));

  // 前下右
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(center.x, boundary_.ymin, boundary_.zmin),
      Position(boundary_.xmax, center.y, center.z)),
    capacity_));

  // 前上左
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(boundary_.xmin, center.y, boundary_.zmin),
      Position(center.x, boundary_.ymax, center.z)),
    capacity_));

  // 前上右
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(center.x, center.y, boundary_.zmin),
      Position(boundary_.xmax, boundary_.ymax, center.z)),
    capacity_));

  // 后下左
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(boundary_.xmin, boundary_.ymin, center.z),
      Position(center.x, center.y, boundary_.zmax)),
    capacity_));

  // 后下右
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(center.x, boundary_.ymin, center.z),
      Position(boundary_.xmax, center.y, boundary_.zmax)),
    capacity_));

  // 后上左
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(Position(boundary_.xmin, center.y, center.z),
      Position(center.x, boundary_.ymax, boundary_.zmax)),
    capacity_));

  // 后上右
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(center, boundary_.max()), capacity_));

  divided_ = true;
}
bool OctreeNode::rayAABBIntersect(const Position& origin,
  const Position& direction, const BoundingBox& aabb) const
{
  Position invDir =
    Position(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z);

  double t1 = (aabb.xmin - origin.x) * invDir.x;
  double t2 = (aabb.xmax - origin.x) * invDir.x;
  double t3 = (aabb.ymin - origin.y) * invDir.y;
  double t4 = (aabb.ymax - origin.y) * invDir.y;
  double t5 = (aabb.zmin - origin.z) * invDir.z;
  double t6 = (aabb.zmax - origin.z) * invDir.z;

  double tmin =
    std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
  double tmax =
    std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

  // 如果tmax < 0，射线在AABB后面
  if (tmax < 0) {
    return false;
  }

  // 如果tmin > tmax，射线不与AABB相交
  if (tmin > tmax) {
    return false;
  }

  return true;
}
void OctreeNode::printTree(int depth) const
{
  // 创建缩进字符串
  std::string indent(depth * 2, ' ');

  // 打印当前节点信息
  std::cout << indent << "└─ Node [depth=" << depth
            << ", spheres=" << spheres_indexs_.size()
            << ", divided=" << (divided_ ? "true" : "false") << "]\n";

  // 打印边界框信息
  Position center = boundary_.getCenter();
  Position size = boundary_.getSize();
  std::cout << indent << "   Bounds: min(" << boundary_.xmin << ", "
            << boundary_.ymin << ", " << boundary_.zmin << ")"
            << " max(" << boundary_.xmax << ", " << boundary_.ymax << ", "
            << boundary_.zmax << ")\n";
  std::cout << indent << "   Center: (" << center.x << ", " << center.y << ", "
            << center.z << ")"
            << " Size: (" << size.x << ", " << size.y << ", " << size.z
            << ")\n";

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