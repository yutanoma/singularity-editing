#include "./convert_zero_face.h"

namespace rp {
  void convert_zero_face(const Eigen::MatrixX3d &_V, const Eigen::MatrixX3i &_F,
                         Eigen::MatrixX3d &V, Eigen::MatrixX3i &F) {
    V = _V;
    F = _F;

    std::vector<std::vector<int>> boundaryLoops;
    igl::boundary_loop(_F, boundaryLoops);

    Eigen::VectorXi isBoundary(_V.rows());
    isBoundary.setZero();

    for (int i = 0; i < boundaryLoops.size(); i++) {
      for (int j = 0; j < boundaryLoops[i].size(); j++) {
        isBoundary(boundaryLoops[i][j]) = 1;
      }
    }

    for (int i = 0; i < _F.rows(); i++) {
      int boundaryCount = 0;
      for (int j = 0; j < 3; j++) {
        if (isBoundary(F(i, j))) {
          boundaryCount++;
        }
      }

      if (boundaryCount >= 2) {
        // このfaceとf0を入れ替える
        auto f0row = F.row(0);
        for (int j = 0; j < 3; j++) {
          F(0, j) = F(i, j);
          F(i, j) = f0row(j);
        }
        return;
      }
    }
  }
}
