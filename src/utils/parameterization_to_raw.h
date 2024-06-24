#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <igl/PI.h>

namespace rp {
  void parameterization_to_raw(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                               const Eigen::MatrixXi &EV,
                               const std::vector<std::vector<int>> &VE,
                               const std::vector<std::vector<double>> &VEAngles,
                               const Eigen::VectorXd &vectorAnglesPerVertex,
                               Eigen::MatrixX3d &vectorField);
}
