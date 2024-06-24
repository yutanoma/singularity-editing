#pragma once

#include <Eigen/Core>
#include "./vector_utils.h"

namespace rp {
  void face_angles(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &normals, Eigen::MatrixX3d &faceAngles);
}
