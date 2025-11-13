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
public:
  BoundingBox boundary_;                // 节点边界
  double min_size_;                     // 最小分割尺寸
  int capacity_;                        // 节点容量
  std::vector<int32_t> spheres_indexs_; // 存储的球体的token（带正负号）
  std::vector<std::unique_ptr<OctreeNode>> children_; // 子节点
  bool divided_;                                      // 是否已分割
  uint64_t morton_code_;                              // Morton编码
  int depth_;                                         // 节点深度
  mutable std::array<std::vector<OctreeNode*>, 6>
    neighborsByFace; // 按面存储的邻居节点指针

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
  std::pair<int32_t, double> queryRay_morton_code(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // 邻居列表版本
  std::pair<int32_t, double> queryRay_Neighbor_search(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // 节点定位版本（无邻域搜索加速）
  std::pair<int32_t, double> queryRay_leaf_find(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // 旧版本的queryRay函数，保留以备对比(后续可删除)
  std::pair<int32_t, double> queryRayold(const Position& origin,
    const Position& direction, int32_t on_surface) const;

  // 打印八叉树结构（用于调试）
  void printTree(int depth = 0, bool showAll = false) const;

  // 新增：查找包含指定球体的所有叶子节点（用于调试）
  std::vector<const OctreeNode*> findLeafNodesContainingSphere(
    int32_t sphere_token) const;
  // 获取Morton编码
  uint64_t getMortonCode() const { return morton_code_; }

  // 获取节点深度
  int getDepth() const { return depth_; }

  // 新增：查询位置所在的叶子节点
  const OctreeNode* findLeafNode(const Position& point) const;

  // 新增：计算射线与当前节点边界的出口距离
  std::pair<BoxFace, double> getExitDistance(
    const Position& origin, const Position& direction) const;
  // 新增：根据出口面推导下一个节点的Morton编码以及位深度
  std::pair<uint64_t, int> deriveNextNodeMorton(BoxFace exit_face) const;

  // 获取特定方向的邻居列表
  std::vector<OctreeNode*>& getNeighbors(BoxFace face)
  {
    return neighborsByFace[static_cast<size_t>(face)];
  }

  // 获取特定方向的邻居列表 (const 版本)
  const std::vector<OctreeNode*>& getNeighbors(BoxFace face) const
  {
    return neighborsByFace[static_cast<size_t>(face)];
  }

  // 添加一个邻居到指定方向
  void addNeighbor(BoxFace face, OctreeNode* neighbor)
  {
    getNeighbors(face).push_back(neighbor);
  }

private:
  // 分割节点
  void subdivide();
  bool shouldSubdivide() const;
  // 新增：Morton编码相关方法
  uint64_t computeChildMortonCode(
    int child_index) const; // 计算子节点的Morton编码
};
void buildLeafMap(OctreeNode* node);
namespace model {
extern std::unordered_map<uint64_t, OctreeNode*> leaf_nodes_map;
} // namespace model

} // namespace openmc
#endif // OCTREE_H
