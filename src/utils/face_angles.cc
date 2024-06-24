#include "./face_angles.h"

namespace rp {
void face_angles(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &normals, Eigen::MatrixX3d &faceAngles) {
  faceAngles.resize(F.rows(), F.cols());
  for (int i = 0; i < F.rows(); i++) {
    Eigen::Vector3d normal = normals.row(i).transpose();
    for (int j = 0; j < F.cols(); j++) {
      int vid = F(i, j);
      int nextVid = F(i, (j + 1) % 3);
      int prevVid = F(i, (j + 2) % 3);
      Eigen::Vector3d base = (V.row(nextVid) - V.row(vid)).transpose();
      Eigen::Vector3d angled = (V.row(prevVid) - V.row(vid)).transpose();

      faceAngles(i, j) = std::abs(rp::get_angle(base, angled, normal));
    }
  }
}
}
