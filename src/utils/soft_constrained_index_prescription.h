#pragma once

#include "index_prescription.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <igl/PI.h>

namespace rp {
  void soft_constrained_index_prescription(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                           const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                                           const Eigen::VectorXi &innerEdges,
                                           const Eigen::VectorXd &angleJumpsPerEdge,
                                           const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                                           const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                                           const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                                           std::vector<Eigen::Triplet<double>> &hardConstraintMatrixTriplets,
                                           const std::vector<double> &hardConstraintVals,
                                           std::vector<Eigen::Triplet<double>> &softConstraintMatrixTriplets,
                                           const std::vector<double> &softConstraintVals,
                                           Eigen::SparseMatrix<double> &matrix,
                                           Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver,
                                           Eigen::VectorXd &cycleCurvature,
                                           Eigen::VectorXd &cycleIndices);
}
