#include "rotation_to_scalar.h"

#include "vector_utils.h"

#include <iostream>

namespace rp {
  void rotation_to_scalar(const Eigen::MatrixX3d &V, const Eigen::MatrixXi &EV, const Eigen::VectorXd &oneForm,
                          const double &globalRotation, const Eigen::VectorXi &treeFathers,
                          const bool &roundPi,
                          Eigen::SparseMatrix<double> &zeroForm2TreeEdges,
                          Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &ldltSolver,
                          bool &needsFactorization, Eigen::VectorXd &periodicAngle) {
    Eigen::VectorXd innerOneForm(V.rows() - 1);
    innerOneForm.setConstant(0);
    int n = 0;
    for (int i = 0; i < V.rows(); i++) {
      if (treeFathers(i) != -1 && n < V.rows() - 1) {
        innerOneForm(n) = oneForm(treeFathers(i));
        n++;
      }
    }

    assert(n == V.rows() - 1);

    // 1-formはEV(i, 1)-EV(i, 0)の値が入っている
    if (needsFactorization || zeroForm2TreeEdges.rows() != innerOneForm.rows() || zeroForm2TreeEdges.cols() != V.rows()) {
      zeroForm2TreeEdges = Eigen::SparseMatrix<double>(innerOneForm.rows(), V.rows());

      // 木の中にあるedgeの数はV.rows() - 1。
      // 各ノードが一つ親を持つため

      std::vector<Eigen::Triplet<double>> triplets = {};

      int nums = 0;
      for (int i = 0; i < treeFathers.rows(); i++) {
        if (treeFathers(i) != -1 && nums < V.rows() - 1) {
          int edgeId = treeFathers(i);

          triplets.emplace_back(Eigen::Triplet<double>(nums, EV(edgeId, 1), 1));
          triplets.emplace_back(Eigen::Triplet<double>(nums, EV(edgeId, 0), -1));
          nums++;
        }
      }

      assert(nums == V.rows() - 1);

      std::cout << "l58" << std::endl;

      zeroForm2TreeEdges.setFromTriplets(triplets.begin(), triplets.end());
      ldltSolver.compute(zeroForm2TreeEdges.transpose() * zeroForm2TreeEdges);

      needsFactorization = false;
    }

    std::cout << "l53" << std::endl;

    periodicAngle = ldltSolver.solve(zeroForm2TreeEdges.transpose() * innerOneForm);

    double diff = globalRotation - periodicAngle(0);

    if (roundPi) {
      for (int i = 0; i < periodicAngle.rows(); i++) {
        periodicAngle(i) = rp::round_pi(periodicAngle(i) + diff);
      }
    }
  }
}
