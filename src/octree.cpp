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
  // Check if sphere is within node boundary
  if (!boundary_.intersects(sphere_token)) {
    return false;
  }

  // If node is not full, insert directly
  if (!divided_ && spheres_indexs_.size() < capacity_) {
    spheres_indexs_.push_back(sphere_token);
    return true;
  }

  // If node is full and not divided, subdivide first
  if (!divided_) {
    // If size is already smaller than minimum subdivision size, don't subdivide, insert directly into current node
    if (!shouldSubdivide()) {
      spheres_indexs_.push_back(sphere_token);
      return true;
    }
    subdivide();

    // Important: redistribute current node's spheres to child nodes
    auto old_spheres = std::move(spheres_indexs_);
    spheres_indexs_.clear();

    for (const auto& token : old_spheres) {
      // Try to insert into child nodes
      bool inserted = false;
      for (auto& child : children_) {
        if (child->insert(token)) {
          inserted = true;
        }
      }
      // If sphere cannot be inserted into any child node, report error
      if (!inserted) {
        fatal_error("Error in OctreeNode::insert: could not reinsert existing "
                    "sphere {} into child nodes.",
          model::surfaces[abs(token) - 1]->id_);
      }
    }
  }

  // Try to insert new sphere into child nodes
  bool inserted_to_child = false;
  for (auto& child : children_) {
    if (child->insert(sphere_token)) {
      inserted_to_child = true;
    }
  }

  // If sphere cannot be inserted into any child node, keep it in current node
  if (!inserted_to_child) {
    fatal_error("Error in OctreeNode::insert: could not reinsert existing "
                "sphere {} into child nodes.",
      model::surfaces[abs(sphere_token) - 1]->id_);
  }

  return true;
}

int32_t OctreeNode::queryPoint(const Position& point) const
{
  // If point is not within node boundary, return -1
  if (!boundary_.contains(point)) {
    return -1;
  }

  // Check spheres in current node
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

  // Check child nodes
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
  Position current_position = origin; // Step (1): set current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  uint64_t old_morton_code = -1;
  bool if_stuck = false;
  int32_t stuck_count = 0;
  BoxFace exit_face = BoxFace::NONE;
  double total_distance = 0.0;
  const OctreeNode* leaf_node = nullptr;

  while (true) {
    bool found = false; // Flag indicating whether it was found
    // First search using Morton code
    if (old_morton_code != -1 && leaf_node != nullptr) {
      auto [next_morton_code, next_depth] =
        leaf_node->deriveNextNodeMorton(exit_face);

      // Loop back through depths until finding corresponding leaf node or cannot go back
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
      // If not found, use position search
      // Find leaf node containing current position
      if (leaf_node == nullptr) {
        leaf_node = findLeafNode(current_position);
      } else {
        leaf_node = leaf_node->findLeafNode(current_position);
      }
      if (leaf_node == nullptr ||
          leaf_node->getMortonCode() == old_morton_code) {
        // If not in any leaf node, move forward a small distance
        current_position += direction * FP_COINCIDENT;
        // Search globally again, if still not found then exit
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
    // If not stuck in the same leaf node, search for spheres in the leaf node containing current position
    if (!if_stuck) {
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // Record sphere token
        }
      }
    }
    // Calculate exit distance of ray with current node boundary
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // If exit distance is infinite, ray will not enter other nodes, exit loop
    exit_face = exit_face1;
    if (exit_distance == INFTY) {
      break;
    }
    // Step (5): update current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    total_distance += exit_distance + FP_COINCIDENT;
    // If exceeds maximum distance, exit loop
    if (total_distance > max_distance) {
      break;
    }
    // Update old_morton_code
    old_morton_code = leaf_node->getMortonCode();
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRay_Neighbor_search(
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // Step (1): set current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  double total_distance = 0.0;
  BoxFace exit_face = BoxFace::NONE;

  const OctreeNode* leaf_node = nullptr; // Leaf node containing current position

  while (true) {
    bool found = false;
    bool if_stuck = false;
    // If exit is not empty and old node exists
    if (exit_face != BoxFace::NONE && old_leaf_node != nullptr) {
      auto& Neighbor_nodes = old_leaf_node->getNeighbors(exit_face);
      // If neighbor list is not empty
      if (!Neighbor_nodes.empty()) {
        // Check if any neighbor node contains current_position
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
      // Find leaf node containing current position
      leaf_node = findLeafNode(current_position);
      if (leaf_node == nullptr || leaf_node == old_leaf_node) {
        // If not in any leaf node, move forward a small distance
        current_position += direction * FP_COINCIDENT;
        // Search again, if still not found then exit
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
        // If not stuck, add new node to old node's neighbor list
        if (old_leaf_node != nullptr) {
          old_leaf_node->addNeighbor(
            exit_face, const_cast<OctreeNode*>(leaf_node));
        }
      }
    }

    // If not stuck in the same leaf node, search for spheres in the leaf node containing current position
    if (!if_stuck) {
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // Record sphere token
        }
      }
    }
    // Calculate exit distance of ray with current node boundary
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    exit_face = exit_face1;
    // If exit distance is infinite, ray will not enter other nodes, exit loop

    if (exit_distance == INFTY) {
      break;
    }
    // Step (5): update current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    total_distance += exit_distance + FP_COINCIDENT;
    // If exceeds maximum distance, exit loop
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
  Position current_position = origin; // Step (1): set current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  bool if_stuck = false;
  int32_t stuck_count = 0;

  while (true) {
    // Find leaf node containing current position
    const OctreeNode* leaf_node = findLeafNode(current_position);
    if (leaf_node == nullptr || leaf_node == old_leaf_node) {
      // If not in any leaf node, move forward a small distance
      current_position += direction * FP_COINCIDENT;
      // Search again, if still not found then exit
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
    // If not stuck in the same leaf node, search for spheres in the leaf node containing current position
    if (!if_stuck) {
      stuck_count = 0;
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // Record sphere token
        }
      }
    }
    // Calculate exit distance of ray with current node boundary
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // If exit distance is infinite, ray will not enter other nodes, exit loop
    if (exit_distance == INFTY) {
      break;
    }
    // Step (5): update current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    old_leaf_node = const_cast<OctreeNode*>(leaf_node);
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRay_leaf_find_old(
  // Without location acceleration
  const Position& origin, const Position& direction, int32_t on_surface,
  double max_distance) const
{
  Position current_position = origin; // Step (1): set current_position=origin
  double minT = INFTY;
  int32_t result_sphere = std::numeric_limits<int32_t>::max();
  OctreeNode* old_leaf_node = nullptr;
  bool if_stuck = false;
  int32_t stuck_count = 0;

  while (true) {
    // Find leaf node containing current position
    const OctreeNode* leaf_node = findLeafNode_old(current_position);
    if (leaf_node == nullptr || leaf_node == old_leaf_node) {
      // If not in any leaf node, move forward a small distance
      current_position += direction * FP_COINCIDENT;
      // Search again, if still not found then exit
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
    // If not stuck in the same leaf node, search for spheres in the leaf node containing current position
    if (!if_stuck) {
      stuck_count = 0;
      for (const auto& sphere_token : leaf_node->spheres_indexs_) {
        bool coincident {std::abs(sphere_token) == std::abs(on_surface)};
        double t = model::surfaces[abs(sphere_token) - 1]->distance(
          origin, direction, coincident);
        if (t > 0 && t < minT) {
          minT = t;
          result_sphere = sphere_token; // Record sphere token
        }
      }
    }
    // Calculate exit distance of ray with current node boundary
    auto [exit_face1, exit_distance] =
      leaf_node->getExitDistance(current_position, direction);
    // If exit distance is infinite, ray will not enter other nodes, exit loop
    if (exit_distance == INFTY) {
      break;
    }
    // Step (5): update current_position
    current_position += direction * (exit_distance + FP_COINCIDENT);
    old_leaf_node = const_cast<OctreeNode*>(leaf_node);
  }
  return {result_sphere, minT};
}

std::pair<int32_t, double> OctreeNode::queryRayold(
  const Position& origin, const Position& direction, int32_t on_surface) const
{
  // If ray does not intersect node boundary, return invalid result
  if (!boundary_.rayIntersect(origin, direction)) {
    return {
      std::numeric_limits<int32_t>::max(), std::numeric_limits<double>::max()};
  }

  std::pair<int32_t, double> result = {
    std::numeric_limits<int32_t>::max(), std::numeric_limits<double>::max()};
  double minT = std::numeric_limits<double>::max();

  // Check spheres in current node
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

  // Check child nodes
  if (divided_) {
    for (const auto& child : children_) {
      auto childResult = child->queryRayold(origin, direction, on_surface);
      if (childResult.first != std::numeric_limits<int32_t>::max() &&
          childResult.second < minT) {
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

  // Uniformly use min/max/center to define boundaries, ensuring no overlap or omission
  // Calculate Morton code for each child node, right-handed coordinate system
  // Positive direction is 1, negative direction is 0, encoding order is zyx

  // Index 0: left-bottom-front (000)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, min.z),
                                   Position(center.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(0), child_depth));

  // Index 1: right-bottom-front (001)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(center.x, min.y, min.z),
                                   Position(max.x, center.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(1), child_depth));

  // Index 2: left-top-front (010)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, center.y, min.z),
                                   Position(center.x, max.y, center.z)),
      min_size_, capacity_, computeChildMortonCode(2), child_depth));

  // Index 3: right-top-front (011)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, min.z), Position(max.x, max.y, center.z)),
    min_size_, capacity_, computeChildMortonCode(3), child_depth));

  // Index 4: left-bottom-back (100)
  children_.push_back(
    std::make_unique<OctreeNode>(BoundingBox(Position(min.x, min.y, center.z),
                                   Position(center.x, center.y, max.z)),
      min_size_, capacity_, computeChildMortonCode(4), child_depth));

  // Index 5: right-bottom-back (101)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, min.y, center.z), Position(max.x, center.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(5), child_depth));

  // Index 6: left-top-back (110)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(min.x, center.y, center.z), Position(center.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(6), child_depth));

  // Index 7: right-top-back (111)
  children_.push_back(std::make_unique<OctreeNode>(
    BoundingBox(
      Position(center.x, center.y, center.z), Position(max.x, max.y, max.z)),
    min_size_, capacity_, computeChildMortonCode(7), child_depth));

  divided_ = true;
}

void OctreeNode::printBasicInfo() const
{
  // Print maximum depth and leaf node count of current tree structure
  int max_depth = 0;
  int total_leaves = 0;
  int total_particles = 0;

  std::function<void(const OctreeNode*, int)> traverse;
  traverse = [&](const OctreeNode* node, int depth) {
    if (!node->divided_) {
      total_leaves++;
      total_particles += node->spheres_indexs_.size();
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

  // Calculate average
  double avg_particles = total_leaves > 0
                           ? static_cast<double>(total_particles) / total_leaves
                           : 0.0;

  // Output basic information
  write_message(fmt::format(
    "Octree basic info: max depth={}, max capacity={}, total leaves={}, total "
    "particles={}, avg particles per leaf={:.2f}",
    max_depth, capacity_, total_leaves, total_particles, avg_particles));
}

void OctreeNode::printTree(int depth, bool showAll) const
{
  if (!showAll && spheres_indexs_.empty() && !divided_) {
    return; // If node is empty and not divided, and not required to show all nodes, skip
  }
  // Create indent string
  std::string indent(depth * 2, ' ');
  // Print current node information
  std::cout << indent << "└─ Node [depth=" << depth
            << ", spheres=" << spheres_indexs_.size()
            << ", divided=" << (divided_ ? "true" : "false") << "]\n";
  if (showAll) {

    // Print bounding box information
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

  // Print spheres in current node
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

  // Recursively print child nodes
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

  // Check if current node contains this sphere
  bool contains_sphere = false;
  for (const auto& token : spheres_indexs_) {
    if (token == sphere_token) {
      contains_sphere = true;
      break;
    }
  }

  // If current node is a leaf node and contains this sphere, add to result
  if (!divided_ && contains_sphere) {
    result.push_back(this);
  }

  // If node is divided, recursively check child nodes
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
  // If point is not within node boundary, return nullptr
  if (!boundary_.contains(point)) {
    return nullptr;
  }

  const OctreeNode* current = this;

  // Iterate downward to find leaf node
  while (current->divided_) {
    Position center = current->boundary_.getCenter();

    // Determine child node index based on point position relative to center
    int child_index = 0;
    if (point.x >= center.x)
      child_index |= 1; // x-axis: 0=left, 1=right
    if (point.y >= center.y)
      child_index |= 2; // y-axis: 0=bottom, 1=top
    if (point.z >= center.z)
      child_index |= 4; // z-axis: 0=front, 1=back

    // Directly access corresponding child node
    current = current->children_[child_index].get();

    // Safety check: if child node is null, return current node
    if (current == nullptr) {
      return const_cast<OctreeNode*>(this);
    }
  }

  return current;
}

const OctreeNode* OctreeNode::findLeafNode_old(const Position& point) const
{
  // If point is not within node boundary, return nullptr
  if (!boundary_.contains(point)) {
    return nullptr;
  }

  const OctreeNode* current = this;

  // Iterate downward to find leaf node
  while (current->divided_) {
    // Brute force: check each child node one by one
    bool found = false;

    // Check all 8 child nodes to find the one containing this point
    for (int i = 0; i < 8; i++) {
      if (current->children_[i] != nullptr &&
          current->children_[i]->boundary_.contains(point)) {
        current = current->children_[i].get();
        found = true;
        break;
      }
    }

    // If no child node containing this point is found, return current node
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
  // Child node's Morton code = parent node code left-shifted by 3 bits + child node index
  return (morton_code_ << 3) | (child_index & 0x7);
}

bool OctreeNode::shouldSubdivide() const
{
  // Check if boundary size is greater than minimum subdivision size
  Position size = boundary_.getSize();
  return (size.x > min_size_ && size.y > min_size_ && size.z > min_size_);
}

std::pair<uint64_t, int> OctreeNode::deriveNextNodeMorton(
  BoxFace exit_face) const
{
  uint64_t next_morton_code = morton_code_;

  // Adjust Morton code directly based on exit face
  switch (exit_face) {
  case BoxFace::MIN_X:
    // X coordinate minus 1: need to find the lowest x bit that is 1, set it to 0, and set all lower x bits to 1
    for (int i = 0; i < depth_; ++i) {
      uint64_t x_bit_pos = 3 * i;
      if ((next_morton_code >> x_bit_pos) & 0x1) {
        // Find first x bit that is 1, set it to 0
        next_morton_code &= ~(1ULL << x_bit_pos);
        // Set all lower x bits to 1
        for (int j = 0; j < i; ++j) {
          next_morton_code |= (1ULL << (3 * j));
        }
        break;
      }
    }
    break;

  case BoxFace::MAX_X:
    // X coordinate plus 1: need to find the lowest x bit that is 0, set it to 1, and set all lower x bits to 0
    for (int i = 0; i < depth_; ++i) {
      uint64_t x_bit_pos = 3 * i;
      if (!((next_morton_code >> x_bit_pos) & 0x1)) {
        // Find first x bit that is 0, set it to 1
        next_morton_code |= (1ULL << x_bit_pos);
        // Set all lower x bits to 0
        for (int j = 0; j < i; ++j) {
          next_morton_code &= ~(1ULL << (3 * j));
        }
        break;
      }
    }
    break;

  case BoxFace::MIN_Y:
    // Y coordinate minus 1: operate on y bit (position 3*i+1)
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
    // Y coordinate plus 1: operate on y bit (position 3*i+1)
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
    // Z coordinate minus 1: operate on z bit (position 3*i+2)
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
    // Z coordinate plus 1: operate on z bit (position 3*i+2)
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