#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "./vector_utils.h"

namespace rp {
  void hodge_star(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV, const Eigen::VectorXi &innerEdges, Eigen::SparseMatrix<double> &hodgeStar);
}