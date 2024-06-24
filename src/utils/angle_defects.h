#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

namespace rp {
  void angle_defects(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                     const Eigen::VectorXd &anglesSumPerVertex,
                     const Eigen::MatrixX3d &faceAngles,
                     const Eigen::SparseMatrix<double> &basisLoops,
                     const std::vector<std::vector<int>> &boundaryLoops,
                     const Eigen::VectorXd &angleJumps,
                     Eigen::VectorXd &angleDefects);
}
