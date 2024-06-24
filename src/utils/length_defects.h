#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
void length_defects(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                    const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                    const Eigen::MatrixXi &EV, const std::vector<std::vector<int>> &VE,
                    const std::vector<std::vector<double>> &VEAngles,
                    const Eigen::VectorXd &vectorAnglesPerVertex,
                    const Eigen::VectorXd &rescalings,
                    const double &velocity,
                    Eigen::VectorXd &projection);
}
