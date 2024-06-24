#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

namespace rp {
void rotation_to_parameterization(const Eigen::MatrixX3d &V, const Eigen::MatrixXi &EV, const Eigen::VectorXd &oneForm,
                                  const double &globalRotation,
                                  Eigen::VectorXd &periodicAngle);
}
