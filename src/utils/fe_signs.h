#pragma once

#include <Eigen/Core>

namespace rp {
  void fe_signs(Eigen::MatrixX3i &F, const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EV, Eigen::MatrixXi &FESigns);
}
