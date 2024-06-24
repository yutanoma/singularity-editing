#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
namespace vertex_edge_angles {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi EV,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
             std::vector<std::vector<int>> &VE, std::vector<std::vector<double>> &VEAngles);

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi EV,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
             std::vector<std::vector<int>> &VE, 
             // VEAnglesは2πに丸められてる
             std::vector<std::vector<double>> &VEAngles, Eigen::VectorXd &anglesSumPerVertex);
}
}
