#pragma once

#include <Eigen/Core>
#include <vector>

// 0-formのparameterizationから、特定の値の等高線を取り出す

namespace rp {
namespace level_set_boundary {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &parameterization,
             const double &value, std::vector<Eigen::MatrixX3d> &paths);
}
}