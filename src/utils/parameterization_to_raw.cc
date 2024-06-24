#include "./parameterization_to_raw.h"

#include <iostream>

namespace rp {
void parameterization_to_raw(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                             const Eigen::MatrixXi &EV,
                             const std::vector<std::vector<int>> &VE,
                             const std::vector<std::vector<double>> &VEAngles,
                             const Eigen::VectorXd &vectorAnglesPerVertex,
                             Eigen::MatrixX3d &vectorField) {
  vectorField.resize(V.rows(), 3);
  vectorField.setZero();
  for (int i = 0; i < vectorAnglesPerVertex.rows(); i++) {
    int vid = i;
    double angle = vectorAnglesPerVertex(i);
    if (angle < 0) {
      angle += 2 * igl::PI;
    }

    for (int j = 0; j < VE[i].size(); j++) {
      double angle0 = VEAngles[i][j], angle1 = VEAngles[i][(j + 1) % VEAngles[i].size()];
      if (angle0 < 0) {
        angle0 += 2 * igl::PI;
      }
      if (angle1 < 0 || j == VE[i].size() - 1) {
        // VEAngles[0]は必ず0なので2πにする
        angle1 += 2 * igl::PI;
      }

      assert(angle0 < angle1);

      if (angle0 <= angle && angle <= angle1) {
        int e0id = VE[i][j];
        int e1id = VE[i][(j + 1) % VEAngles[i].size()];
        int v0id = EV(e0id, 0) == vid ? EV(e0id, 1) : EV(e0id, 0);
        int v1id = EV(e1id, 0) == vid ? EV(e1id, 1) : EV(e1id, 0);

        Eigen::Vector3d e0Vec = (V.row(v0id) - V.row(vid)).transpose().normalized();
        Eigen::Vector3d e1Vec = (V.row(v1id) - V.row(vid)).transpose().normalized();
        Eigen::Vector3d normal = e0Vec.cross(e1Vec).normalized();

        double rotateAngle = angle - angle0;

        vectorField.row(i) = Eigen::AngleAxisd(rotateAngle, normal) * e0Vec;
        break;
      }
      if (VE[i].size() - 1 == j)
        std::cout << i << ", " << angle << ", " << angle0 << ", " << angle1 << std::endl;

      assert(j != VE[i].size() - 1);
    }
  }
}
}
