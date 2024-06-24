#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace rp {
namespace non_dual_cycles {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE,
             const Eigen::MatrixXi &FE,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             Eigen::SparseMatrix<double>& basisCycles);

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
              const std::vector<std::vector<int>> &VE,
              const Eigen::MatrixXi &FE,
              const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
              const std::vector<std::vector<int>> &boundaryLoops,
              Eigen::SparseMatrix<double> & basisCycles,
              std::vector<Eigen::Triplet<double>> &basisCycleTriplets);
}
}
