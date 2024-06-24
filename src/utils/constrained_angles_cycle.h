#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>

#include "./edge_paths.h"

namespace rp {
void constrained_angles_cycle(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                              const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                              const std::vector<std::vector<int>> &VE,
                              const std::vector<std::vector<double>> &VEAngles,
                              const Eigen::VectorXi &treeFathers,
                              const Eigen::VectorXd &angleJumps,
                              const std::vector<std::vector<int>> &boundaryLoops,
                              const std::vector<std::vector<int>> &outflowPathVids,
                              const std::vector<std::vector<int>> &inflowPathVids,
                              const std::vector<std::vector<EdgePath>> &brushPaths,
                              const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                              const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                              const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                              // ベクトル場とparameterizationの角度差
                              const double &rotation,
                              double &globalRotation,
                              Eigen::VectorXd &angleDefects,
                              Eigen::VectorXd &cycleIndices,
                              Eigen::SparseMatrix<double> &matrix,
                              Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver);
}
