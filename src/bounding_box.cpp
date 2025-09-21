#include <algorithm> // for min, max

#include "openmc/bounding_box.h"
#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/surface.h"

namespace openmc {

bool BoundingBox::contains(const Position& point) const
{
  return point.x >= xmin && point.x <= xmax && point.y >= ymin &&
         point.y <= ymax && point.z >= zmin && point.z <= zmax;
}

bool BoundingBox::intersects(const int32_t& sphere_token) const
{
  vector<double> triso_center =
    model::surfaces[abs(sphere_token) - 1]->get_center();
  double closestX = std::max(xmin, std::min(triso_center[0], xmax));
  double closestY = std::max(ymin, std::min(triso_center[1], ymax));
  double closestZ = std::max(zmin, std::min(triso_center[2], zmax));

  double distance =
    std::sqrt((closestX - triso_center[0]) * (closestX - triso_center[0]) +
              (closestY - triso_center[1]) * (closestY - triso_center[1]) +
              (closestZ - triso_center[2]) * (closestZ - triso_center[2]));

  return distance <= model::surfaces[abs(sphere_token) - 1]->get_radius();
}

Position BoundingBox::getCenter() const
{
  return Position(
    (xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0);
}

Position BoundingBox::getSize() const
{
  return Position(xmax - xmin, ymax - ymin, zmax - zmin);
}
} // namespace openmc