#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

namespace rp {
  void rotation_to_scalar(const Eigen::MatrixX3d &V, const Eigen::MatrixXi &EV, const Eigen::VectorXd &oneForm,
                          const double &globalRotation, const Eigen::VectorXi &treeFathers,
                          const bool &roundPi,
                          Eigen::SparseMatrix<double> &zeroForm2TreeEdges,
                          Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &ldltSolver,
                          bool &needsFactorization, Eigen::VectorXd &periodicAngle);
}