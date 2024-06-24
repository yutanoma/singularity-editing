#include "./face_based_to_vertex_based.h"

#include <complex>
#include <iostream>

namespace rp {
namespace face_based_to_vertex_based {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas, const Eigen::MatrixX3d &B3,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::MatrixX3d &vectorFieldPerFace,
             const Eigen::VectorXd &anglesSum,
             Eigen::VectorXd &vectorAnglesPerVertex) {
  vectorAnglesPerVertex.resize(V.rows());
  vectorAnglesPerVertex.setZero();

  for (int i = 0; i < V.rows(); i++) {
    // 各edgeに対してベクトル場の回転角を求める
    assert(VE[i].size() != 0);

    std::complex<double> sum(.0, .0);

    // if (i < 10)
    //   std::cout << "[" << V.row(i) << "]: ";

    std::vector<int> edgeList = VE[i];
    for (int j = 0; j < edgeList.size(); j++) {
      int edgeId = edgeList[j];

      if (j == edgeList.size() - 1) {
        if (EF(edgeId, 0) == -1 || EF(edgeId, 1) == -1) {
          break;
        }
      }

      int nextEdgeId = edgeList[(j + 1) % edgeList.size()];

      int faceId = -1;

      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) {
          if (EF(edgeId, k) == EF(nextEdgeId, l)) {
            faceId = EF(edgeId, k);
            break;
          }
        }
        if (faceId != -1) {
          break;
        }
      }

      assert(faceId != -1);

      int edgeOtherVid = EV(edgeId, 0) == i ? EV(edgeId, 1) : EV(edgeId, 0);
      Eigen::Vector3d edgeVec = (V.row(edgeOtherVid) - V.row(i)).transpose();
      Eigen::Vector3d fieldVec = vectorFieldPerFace.row(faceId).transpose();
      Eigen::Vector3d normalVec = B3.row(faceId).transpose();
      double _angle = get_angle(edgeVec, fieldVec, normalVec) + VEAngles[i][j];
      double angle = round_pi(_angle * 2 * igl::PI / anglesSum(i));
      if (angle < 0) {
        // [-π, π]を[0, 2π]にする
        angle += 2 * igl::PI;
      }

      double weight = .0;
      for (int k = 0; k < 3; k++) {
        if (F(faceId, k) == i) {
          weight = faceAngles(faceId, k);
        }
      }

      // std::cout << angle << ", " << i << ", " << j << std::endl;

      // if (i < 10) {
      //   // std::cout << edgeList[j] << "(" << V.row(edgeOtherVid) << "), (" << vectorFieldPerFace.row(faceId) << "), " << faceId << ", ";
      //   std::cout << F.row(faceId) << ", " << faceAngles.row(faceId) << ", "
      //             << VEAngles[i][j] << ", " << anglesSum(i) <<  ", " << angle <<  " | ";
      // }

      sum += weight * std::exp(std::complex<double>(.0, angle));
    }

    // if (i < 10)
    //   std::cout << std::endl;

    double arg = std::atan2(sum.imag(), sum.real());
    if (arg < 0) {
      arg += 2 * igl::PI;
    }

    vectorAnglesPerVertex(i) = arg;

    // std::cout << vectorAnglesPerVertex(i) << std::endl;
  }
}
}
}
