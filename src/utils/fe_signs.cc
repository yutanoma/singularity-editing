#include "./fe_signs.h"

namespace rp {
void fe_signs(Eigen::MatrixX3i &F, const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EV, Eigen::MatrixXi &FESigns) {
  FESigns.resize(FE.rows(), FE.cols());
  for (int i = 0; i < FE.rows(); i++) {
    for (int j = 0; j < FE.cols(); j++) {
      int edgeId = FE(i, j);
      for (int k = 0; k < 3; k++) {
        if (EV(edgeId, 0) == F(i, k) && EV(edgeId, 1) == F(i, (k + 1) % 3)) {
          FESigns(i, j) = 1;
          break;
        } else if (EV(edgeId, 1) == F(i, k) && EV(edgeId, 0) == F(i, (k + 1) % 3)) {
          FESigns(i, j) = -1;
          break;
        }
        assert(k != 2);
      }
    }
  }
}
}
