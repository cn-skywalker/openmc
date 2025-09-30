#ifndef OCTREE_H
#define OCTREE_H
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include <openmc/bounding_box.h>
#include <openmc/position.h>
#include <openmc/surface.h>

namespace openmc {

// 八叉树节点类
class OctreeNode {
private:
  BoundingBox boundary_;                // 节点边界
  int capacity_;                        // 节点容量
  std::vector<int32_t> spheres_indexs_; // 存储的球体的token（带正负号）
  std::vector<std::unique_ptr<OctreeNode>> children_; // 子节点
  bool divided_;                                      // 是否已分割

public:
  OctreeNode(const BoundingBox& boundary, int capacity)
    : boundary_(boundary), capacity_(capacity), divided_(false)
  {}
  ~OctreeNode() = default;
  // 插入球体
  bool insert(const int32_t& sphere_token);

  // 查询点所在的球体，返回球体ID（如果不在任何球体内返回-1）
  int32_t queryPoint(const Position& point) const;

  // 查询射线碰到的第一个球体，返回球体ID和交点（如果没有碰到任何球体返回-1）
  std::pair<int, Position> queryRay(const Position& origin,
    const Position& direction, int32_t on_surface) const;

  // 打印八叉树结构（用于调试）
  void printTree(int depth = 0, bool showAll = false) const;

private:
  // 分割节点
  void subdivide();
  // 射线与BoundingBox相交检测
  bool rayAABBIntersect(const Position& origin, const Position& direction,
    const BoundingBox& aabb) const;
};

} // namespace openmc
#endif // OCTREE_H
