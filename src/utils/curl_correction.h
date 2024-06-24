#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
namespace curl_correction {

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixX3d &B1,
             const Eigen::MatrixX3d &B2, const Eigen::MatrixX3d &B3,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::VectorXd &vectorAnglesPerVertex,
             const Eigen::VectorXd &anglesSumPerVertex,
             Eigen::VectorXd &rescalings);
}
}
