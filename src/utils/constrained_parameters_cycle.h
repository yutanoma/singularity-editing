#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>

#include "./edge_paths.h"

namespace rp {
  void constrained_parameters_cycle(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                                    const Eigen::VectorXi &treeFathers,
                                    const Eigen::VectorXd &projectionPerEdge,
                                    const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                                    const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                                    const std::vector<std::vector<EdgePath>> &indexPrescriptionPaths,
                                    const std::vector<double> indexPrescriptionPathIndices,
                                    const std::vector<std::vector<EdgePath>> &brushPaths,
                                    const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                                    const double &levelsetVal,
                                    double &globalRotation,
                                    Eigen::VectorXd &cycleIndices,
                                    Eigen::SparseMatrix<double> &matrix,
                                    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver);
}
