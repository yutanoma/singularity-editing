#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <igl/PI.h>

namespace rp {
  void singularity_decomposition(const Eigen::MatrixXi &EV, const Eigen::VectorXd &projectionPerEdge,
                                 const Eigen::VectorXd &parameterization, const Eigen::SparseMatrix<double> &basisLoopsMatrix,
                                 Eigen::VectorXd &indices);

  void singularity_decomposition(const Eigen::MatrixXi &EV, const Eigen::VectorXd &projectionPerEdge,
                                 const Eigen::VectorXd &parameterization, const Eigen::SparseMatrix<double> &basisLoopsMatrix,
                                 Eigen::VectorXd &oneForm, Eigen::VectorXd &indices);
}
