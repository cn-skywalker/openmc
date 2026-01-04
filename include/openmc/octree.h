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

// Octree ray tracing mode
enum class OctreeRayTraceMode {
  MORTON_CODE,     // Node location based on Morton code
  NEIGHBOR_SEARCH, // Node location based on neighbor list
  LEAF_FIND,       // Node location based on leaf node search
  LEAF_FIND_OLD,   // Old version of leaf node search based location
  QUERY_RAY_OLD    // Old version of queryRay function
};

// Octree node class
class OctreeNode {
public:
  BoundingBox boundary_;                // Node boundary
  double min_size_;                     // Minimum subdivision size
  int capacity_;                        // Node capacity
  std::vector<int32_t> spheres_indexs_; // Stored sphere tokens (with sign)
  std::vector<std::unique_ptr<OctreeNode>> children_; // Child nodes
  bool divided_;                                      // Whether subdivided
  uint64_t morton_code_;                              // Morton code
  int depth_;                                         // Node depth
  mutable std::array<std::vector<OctreeNode*>, 6>
    neighborsByFace; // Neighbor node pointers stored by face
  OctreeNode(const BoundingBox& boundary, double min_size, int capacity,
    uint64_t morton_code = 1, int depth = 0, bool divided = false)
    : boundary_(boundary), min_size_(min_size), capacity_(capacity),
      morton_code_(morton_code), depth_(depth), divided_(divided)
  {}
  ~OctreeNode() = default;
  // Insert sphere
  bool insert(const int32_t& sphere_token);

  // Query sphere containing point, returns sphere ID (returns -1 if not in any sphere)
  int32_t queryPoint(const Position& point) const;

  // Query first sphere hit by ray, returns sphere ID and intersection point (returns -1 if no sphere hit)
  std::pair<int32_t, double> queryRay_morton_code(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // Neighbor list version
  std::pair<int32_t, double> queryRay_Neighbor_search(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // Node location version (without neighbor search acceleration)
  std::pair<int32_t, double> queryRay_leaf_find(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // Old version of queryRay_leaf_find function, kept for comparison (can be deleted later)
  std::pair<int32_t, double> queryRay_leaf_find_old(const Position& origin,
    const Position& direction, int32_t on_surface,
    double max_distance = INFTY) const;

  // Old version of queryRay function, kept for comparison (can be deleted later)
  std::pair<int32_t, double> queryRayold(const Position& origin,
    const Position& direction, int32_t on_surface) const;

  // Print octree structure (for debugging)
  void printTree(int depth = 0, bool showAll = false) const;

  // Print octree basic information
  void printBasicInfo() const;

  // Find all leaf nodes containing specified sphere (for debugging)
  std::vector<const OctreeNode*> findLeafNodesContainingSphere(
    int32_t sphere_token) const;
  // Get Morton code
  uint64_t getMortonCode() const { return morton_code_; }

  // Get node depth
  int getDepth() const { return depth_; }

  // Query leaf node containing position
  const OctreeNode* findLeafNode(const Position& point) const;
  const OctreeNode* findLeafNode_old(const Position& point) const;

  // Calculate exit distance of ray with current node boundary
  std::pair<BoxFace, double> getExitDistance(
    const Position& origin, const Position& direction) const;
  // Derive next node's Morton code and bit depth based on exit face
  std::pair<uint64_t, int> deriveNextNodeMorton(BoxFace exit_face) const;

  // Get neighbor list for specific direction
  std::vector<OctreeNode*>& getNeighbors(BoxFace face)
  {
    return neighborsByFace[static_cast<size_t>(face)];
  }

  // Get neighbor list for specific direction (const version)
  const std::vector<OctreeNode*>& getNeighbors(BoxFace face) const
  {
    return neighborsByFace[static_cast<size_t>(face)];
  }

  // Add a neighbor to specified direction
  void addNeighbor(BoxFace face, OctreeNode* neighbor)
  {
    getNeighbors(face).push_back(neighbor);
  }

private:
  // Subdivide node
  void subdivide();
  bool shouldSubdivide() const;
  // Morton code related methods
  uint64_t computeChildMortonCode(
    int child_index) const; // Calculate child node's Morton code
};
void buildLeafMap(OctreeNode* node);
namespace model {
extern std::unordered_map<uint64_t, OctreeNode*> leaf_nodes_map;
} // namespace model

} // namespace openmc
#endif // OCTREE_H
