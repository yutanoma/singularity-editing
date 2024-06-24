#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
namespace stripe_patterns {

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::VectorXd &vectorAnglesPerVertex,
             const Eigen::VectorXd &rescalings, const double &velocity,
             Eigen::VectorXd &parameterization);
}
}
