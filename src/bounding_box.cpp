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
  double tmin = 0.0; // 射线 t >= 0
  double tmax = std::numeric_limits<double>::max();

  for (int i = 0; i < 3; ++i) {
    double dir = direction[i];
    double minVal = min()[i];
    double maxVal = max()[i];

    if (std::abs(dir) < FP_PRECISION) {
      // 射线平行于该轴
      if (origin[i] < minVal || origin[i] > maxVal) {
        return std::numeric_limits<double>::max(); // 不相交
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
        return std::numeric_limits<double>::max(); // 不相交
      }
    }
  }

  return tmin; // 返回最近的相交距离
}

// 在 bounding_box.cpp 中添加以下实现：

std::pair<BoxFace, double> BoundingBox::rayIntersectionDistances(
  const Position& origin, const Position& direction) const
{
  double tmin = 0.0;
  double tmax = std::numeric_limits<double>::max();
  BoxFace exit_face = BoxFace::NONE; // 初始化出口面

  for (int i = 0; i < 3; ++i) {
    double dir = direction[i];
    double minVal = min()[i];
    double maxVal = max()[i];

    if (std::abs(dir) < FP_PRECISION) {
      // 射线平行于该轴
      if (origin[i] < minVal || origin[i] > maxVal) {
        return {BoxFace::NONE, std::numeric_limits<double>::max()}; // 不相交
      }
    } else {
      double invD = 1.0 / dir;
      double t0 = (minVal - origin[i]) * invD;
      double t1 = (maxVal - origin[i]) * invD;

      if (invD < 0.0) {
        std::swap(t0, t1);
      }

      // 更新 tmin (入口)
      if (t0 > tmin) {
        tmin = t0;
      }

      // 更新 tmax (出口) 并记录是哪个面
      if (t1 < tmax) {
        tmax = t1;
        // 根据当前轴和方向，确定出口面
        if (i == 0) { // X轴
          exit_face = (dir > 0) ? BoxFace::MAX_X : BoxFace::MIN_X;
        } else if (i == 1) { // Y轴
          exit_face = (dir > 0) ? BoxFace::MAX_Y : BoxFace::MIN_Y;
        } else if (i == 2) { // Z轴
          exit_face = (dir > 0) ? BoxFace::MAX_Z : BoxFace::MIN_Z;
        }
      }

      if (tmax <= tmin) {
        return {BoxFace::NONE, std::numeric_limits<double>::max()};
      }
    }
  }

  return {exit_face, tmax};
}

} // namespace openmc