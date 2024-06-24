#include "singularity_decomposition.h"

namespace rp {
void singularity_decomposition(const Eigen::MatrixXi &EV, const Eigen::VectorXd &projectionPerEdge,
                               const Eigen::VectorXd &parameterization, const Eigen::SparseMatrix<double> &basisLoopsMatrix,
                               Eigen::VectorXd &indices) {
  Eigen::VectorXd oneForm;
  singularity_decomposition(EV, projectionPerEdge, parameterization, basisLoopsMatrix, oneForm, indices);
}

void singularity_decomposition(const Eigen::MatrixXi &EV, const Eigen::VectorXd &projectionPerEdge,
                               const Eigen::VectorXd &parameterization, const Eigen::SparseMatrix<double> &basisLoopsMatrix,
                               Eigen::VectorXd &oneForm, Eigen::VectorXd &indices) {
  oneForm.resize(EV.rows());

  for (int i = 0; i < EV.rows(); i++) {
    int v0 = EV(i, 0), v1 = EV(i, 1);
    double diff = parameterization(v1) - parameterization(v0);
    double objective = projectionPerEdge(i);

    if (diff - objective > igl::PI) {
      while (diff - objective > igl::PI) {
        diff -= 2 * igl::PI;
      }
    } else if (objective - diff > igl::PI) {
      while (objective - diff > igl::PI) {
        diff += 2 * igl::PI;
      }
    }

    oneForm(i) = diff;
  }

  Eigen::VectorXd loopSums = basisLoopsMatrix * oneForm;

  indices.resize(EV.rows());

  for (int i = 0; i < loopSums.rows(); i++)
    indices(i) = (int)std::round(loopSums(i) / (2 * igl::PI));
};
}

