#ifndef OPENMC_BOUNDING_BOX_H
#define OPENMC_BOUNDING_BOX_H

#include <algorithm> // for min, max

#include "openmc/constants.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Coordinates for an axis-aligned cuboid that bounds a geometric object.
//==============================================================================

class BoundingBox {
public:
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;

  BoundingBox()
    : xmin(-INFTY), xmax(INFTY), ymin(-INFTY), ymax(INFTY), zmin(-INFTY),
      zmax(INFTY)
  {}

  BoundingBox(const Position& min, const Position& max)
    : xmin(min.x), xmax(max.x), ymin(min.y), ymax(max.y), zmin(min.z),
      zmax(max.z)
  {}

  BoundingBox(double x1, double x2, double y1, double y2, double z1, double z2)
    : xmin(x1), xmax(x2), ymin(y1), ymax(y2), zmin(z1), zmax(z2)
  {}

  inline BoundingBox operator&(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result &= other;
  }

  inline BoundingBox operator|(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  inline BoundingBox& operator&=(const BoundingBox& other)
  {
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  inline BoundingBox& operator|=(const BoundingBox& other)
  {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }

  inline Position min() const { return {xmin, ymin, zmin}; }
  inline Position max() const { return {xmax, ymax, zmax}; }
  bool contains(const Position& point) const;
  bool intersects(const int32_t& sphere_token) const;
  Position getCenter() const;
  Position getSize() const;
  bool rayIntersect(const Position& origin, const Position& direction) const;
  // 计算射线到边界框的最近相交距离
  double rayDistance(const Position& origin, const Position& direction) const;
};

} // namespace openmc

#endif
