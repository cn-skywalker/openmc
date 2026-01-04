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

bool BoundingBox::rayIntersect(
  const Position& origin, const Position& direction) const
{
  return rayDistance(origin, direction) < std::numeric_limits<double>::max();
}

double BoundingBox::rayDistance(
  const Position& origin, const Position& direction) const
{
  double tmin = 0.0; // Ray t >= 0
  double tmax = std::numeric_limits<double>::max();

  for (int i = 0; i < 3; ++i) {
    double dir = direction[i];
    double minVal = min()[i];
    double maxVal = max()[i];

    if (std::abs(dir) < FP_PRECISION) {
      // Ray is parallel to this axis
      if (origin[i] < minVal || origin[i] > maxVal) {
        return std::numeric_limits<double>::max(); // No intersection
      }
      // Otherwise don't restrict t range, skip
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
        return std::numeric_limits<double>::max(); // No intersection
      }
    }
  }

  return tmin; // Return the nearest intersection distance
}

std::pair<BoxFace, double> BoundingBox::rayIntersectionDistances(
  const Position& origin, const Position& direction) const
{
  Position minVal = min();
  Position maxVal = max();

  bool inside = (origin[0] > minVal[0] - FP_PRECISION &&
                 origin[0] < maxVal[0] + FP_PRECISION &&
                 origin[1] > minVal[1] - FP_PRECISION &&
                 origin[1] < maxVal[1] + FP_PRECISION &&
                 origin[2] > minVal[2] - FP_PRECISION &&
                 origin[2] < maxVal[2] + FP_PRECISION);

  if (!inside) {
    // External case: use full algorithm
    double tmin = 0.0;
    double tmax = std::numeric_limits<double>::max();
    BoxFace exit_face = BoxFace::NONE;

    for (int i = 0; i < 3; ++i) {
      double dir = direction[i];

      if (std::abs(dir) < FP_PRECISION) {
        if (origin[i] < minVal[i] - FP_COINCIDENT ||
            origin[i] > maxVal[i] + FP_COINCIDENT) {
          return {BoxFace::NONE, std::numeric_limits<double>::max()};
        }
      } else {
        double invD = 1.0 / dir;
        double t0 = (minVal[i] - origin[i]) * invD;
        double t1 = (maxVal[i] - origin[i]) * invD;

        if (invD < 0.0) {
          std::swap(t0, t1);
        }

        if (t0 > tmin)
          tmin = t0;
        if (t1 < tmax) {
          tmax = t1;
          if (i == 0) {
            exit_face = (dir > 0) ? BoxFace::MAX_X : BoxFace::MIN_X;
          } else if (i == 1) {
            exit_face = (dir > 0) ? BoxFace::MAX_Y : BoxFace::MIN_Y;
          } else {
            exit_face = (dir > 0) ? BoxFace::MAX_Z : BoxFace::MIN_Z;
          }
        }

        if (tmax <= tmin + FP_COINCIDENT) {
          return {BoxFace::NONE, std::numeric_limits<double>::max()};
        }
      }
    }
    return {exit_face, tmax};
  } else {
    // Internal case: improved simplified algorithm
    double tmax = std::numeric_limits<double>::max();
    BoxFace exit_face = BoxFace::NONE;

    for (int i = 0; i < 3; ++i) {
      double dir = direction[i];

      if (std::abs(dir) < FP_PRECISION) {
        continue;
      }

      double invD = 1.0 / dir;
      double t_exit;
      BoxFace candidate_face;

      if (dir > 0) {
        t_exit = (maxVal[i] - origin[i]) * invD;
        if (i == 0)
          candidate_face = BoxFace::MAX_X;
        else if (i == 1)
          candidate_face = BoxFace::MAX_Y;
        else
          candidate_face = BoxFace::MAX_Z;
      } else {
        t_exit = (minVal[i] - origin[i]) * invD;
        if (i == 0)
          candidate_face = BoxFace::MIN_X;
        else if (i == 1)
          candidate_face = BoxFace::MIN_Y;
        else
          candidate_face = BoxFace::MIN_Z;
      }

      // Critical fix: ignore negative values or overly small positive values
      if (t_exit > FP_COINCIDENT && t_exit < tmax) {
        tmax = t_exit;
        exit_face = candidate_face;
      }
    }

    return {exit_face, tmax};
  }
}

} // namespace openmc