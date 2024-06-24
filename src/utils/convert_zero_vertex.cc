#include "./convert_zero_vertex.h"

namespace rp {
  void convert_zero_vertex(const Eigen::MatrixX3d &_V, const Eigen::MatrixX3i &_F,
                         Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, int &originalZeroIdx) {
    V = _V;
    F = _F;

    std::vector<std::vector<int>> boundaryLoops;
    igl::boundary_loop(_F, boundaryLoops);

    Eigen::VectorXi isBoundary(_V.rows());
    isBoundary.setZero();

    for (int i = 0; i < boundaryLoops.size(); i++) {
      for (int j = 3; j < boundaryLoops[i].size(); j++) {
        int originalZeroIdx = boundaryLoops[i][j];
        V.row(0) = _V.row(originalZeroIdx);
        V.row(originalZeroIdx) = _V.row(0);

        for (int k = 0; k < F.rows(); k++) {
          for (int l = 0; l < 3; l++) {
            if (F(k, l) == 0) {
              F(k, l) = originalZeroIdx;
            } else if (F(k, l) == originalZeroIdx) {
              F(k, l) = 0;
            }
          }
        }

        return;
      }
    }
  }
}
