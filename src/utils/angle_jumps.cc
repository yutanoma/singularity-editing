#include "angle_jumps.h"

#include "vector_utils.h"

#include <array>
#include <iostream>

namespace rp {
void angle_jumps(const Eigen::MatrixX3i &F, const Eigen::MatrixXi &EV,
                 const std::vector<std::vector<int>> &VE,
                 const std::vector<std::vector<double>> &VEAngles,
                 Eigen::VectorXd &angleJumpPerEdge) {
  angleJumpPerEdge.resize(EV.rows());
  angleJumpPerEdge.setZero();

  for (int i = 0; i < EV.rows(); i++) {
    int edgeId = i;
    std::array<double, 2> edge2BaseAngles = {0, 0};

    assert(EV(edgeId, 0) < EV(edgeId, 1));

    for (int j = 0; j < 2; j++) {
      int vid = EV(edgeId, j);
      std::vector<int> adjacentEdges = VE[vid];
      double angle = 999999;

      for (int k = 0; k < adjacentEdges.size(); k++) {
        if (adjacentEdges[k] == edgeId) {
          angle = VEAngles[vid][k];
        }
      }

      assert(angle != 999999);

      if (j == 0) {
        edge2BaseAngles[j] = rp::round_pi(igl::PI - angle);
      } else {
        edge2BaseAngles[j] = rp::round_pi(-angle);
      }
    }

    angleJumpPerEdge(i) = edge2BaseAngles[0] - edge2BaseAngles[1];
  }
};
}
