#pragma once

#include <Eigen/Core>
#include <igl/boundary_loop.h>
#include <vector>

namespace rp {
  void convert_zero_face(const Eigen::MatrixX3d &_V, const Eigen::MatrixX3i &_F,
                         Eigen::MatrixX3d &V, Eigen::MatrixX3i &F);
}
