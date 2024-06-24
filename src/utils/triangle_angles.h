#pragma once

#include <igl/PI.h>
#include <Eigen/Core>

namespace rp {
  // 長さがa, b, cの三角形の角度を求める
  void triangle_angles(const double &a, const double &b, const double &c,
                       double &A, double &B, double &C);
}
