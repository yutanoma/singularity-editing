#pragma once

#include <Eigen/Core>

namespace rp {
  void save_developable(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const std::string &filename);
}
