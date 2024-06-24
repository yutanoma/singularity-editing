#pragma once

#include <Eigen/Core>

namespace rp {
  void singularity_spheres(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                           const Eigen::VectorXd &indices, Eigen::MatrixXi &sF,
                           Eigen::MatrixXd &sV, Eigen::MatrixXd &sC, double radius);
}
