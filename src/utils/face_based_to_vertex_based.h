#pragma once

#include "vector_utils.h"

#include <Eigen/Core>
#include <vector>

namespace rp {
namespace face_based_to_vertex_based {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas, const Eigen::MatrixX3d &B3,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::MatrixX3d &vectorFieldPerFace,
             const Eigen::VectorXd &anglesSum,
             Eigen::VectorXd &vectorAnglesPerVertex);
}
}
