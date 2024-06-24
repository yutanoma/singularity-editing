#pragma once

#include <Eigen/Core>

namespace rp {
  void write_ply(const std::string &filename, const Eigen::MatrixX3i &F,
                 const Eigen::MatrixX3d &V, const Eigen::MatrixX3d &C);
}
