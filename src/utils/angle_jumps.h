#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

namespace rp {
  void angle_jumps(const Eigen::MatrixX3i &F, const Eigen::MatrixXi &EV,
                   const std::vector<std::vector<int>> &VE,
                   const std::vector<std::vector<double>> &VEAngles,
                   Eigen::VectorXd &angleJumpPerEdge);
}
