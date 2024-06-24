#include "length_defects.h"

#include <igl/PI.h>

namespace rp {
void length_defects(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                    const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                    const Eigen::MatrixXi &EV, const std::vector<std::vector<int>> &VE,
                    const std::vector<std::vector<double>> &VEAngles,
                    const Eigen::VectorXd &vectorAnglesPerVertex,
                    const Eigen::VectorXd &rescalings,
                    const double &velocity,
                    Eigen::VectorXd &projection) {
  projection.resize(EV.rows());
  for (int i = 0; i < EV.rows(); i++) {
    // EV(i, 0)->EV(i, 1)の方向が正
    int v0 = EV(i, 0), v1 = EV(i, 1);
    Eigen::Vector3d edgeVec = (V.row(v1) - V.row(v0)).transpose();
    double length = edgeVec.norm();

    Eigen::Vector2d theta;
    theta.setZero();
    for (int j = 0; j < 2; j++) {
      int vid = EV(i, j);
      bool visited = false;
      for (int k = 0; k < VE[vid].size(); k++) {
        if (VE[vid][k] == i) {
          theta(j) = VEAngles[vid][k];
          visited = true;
          // v0->v1の向きが正なので、もしvidがv1ならπ回転させる
          if (vid == v1) {
            theta(j) += igl::PI;
          }
          break;
        }
      }
      assert(visited);
    }

    double arg0 = velocity * rescalings(EV(i, 0)) * std::cos(vectorAnglesPerVertex(EV(i, 0)) - theta(0));
    double arg1 = velocity * rescalings(EV(i, 1)) * std::cos(vectorAnglesPerVertex(EV(i, 1)) - theta(1));

    double angle = length * (arg0 + arg1) / 2;

    projection(i) = angle;
  }
}
}
