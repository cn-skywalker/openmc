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
  double min_size_;                     // 最小分割尺寸
  int capacity_;                        // 节点容量
  std::vector<int32_t> spheres_indexs_; // 存储的球体的token（带正负号）
  std::vector<std::unique_ptr<OctreeNode>> children_; // 子节点
  bool divided_;                                      // 是否已分割
  uint64_t morton_code_;                              // Morton编码
  int depth_;                                         // 节点深度

public:
  OctreeNode(const BoundingBox& boundary, double min_size, int capacity,
    uint64_t morton_code = 1, int depth = 0, bool divided = false)
    : boundary_(boundary), min_size_(min_size), capacity_(capacity),
      morton_code_(morton_code), depth_(depth), divided_(divided)
  {}
  ~OctreeNode() = default;
  // 插入球体
  bool insert(const int32_t& sphere_token);

  // 查询点所在的球体，返回球体ID（如果不在任何球体内返回-1）
  int32_t queryPoint(const Position& point) const;

  // 查询射线碰到的第一个球体，返回球体ID和交点（如果没有碰到任何球体返回-1）
  std::pair<int32_t, double> queryRay(const Position& origin,
    const Position& direction, int32_t on_surface) const;

  // 打印八叉树结构（用于调试）
  void printTree(int depth = 0, bool showAll = false) const;

  // 获取Morton编码
  uint64_t getMortonCode() const { return morton_code_; }

  // 获取节点深度
  int getDepth() const { return depth_; }

private:
  // 分割节点
  void subdivide();
  bool shouldSubdivide() const;
  // 新增：Morton编码相关方法
  uint64_t computeChildMortonCode(
    int child_index) const; // 计算子节点的Morton编码
  static uint64_t computeMortonCode(const Position& pos, const Position& min,
    const Position& max, int max_depth); // 根据位置计算Morton编码
  static void decodeMortonCode(uint64_t code, int depth, Position& min,
    Position& max, const Position& root_min,
    const Position& root_max); // 从Morton编码解码位置范围
};

} // namespace openmc
#endif // OCTREE_H
