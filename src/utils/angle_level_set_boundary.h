#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

namespace rp {
namespace angle_level_set_boundary {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &zeroForm,
             const Eigen::VectorXd &oneForm, const Eigen::MatrixXi &FESigns,
             const double &value, std::vector<Eigen::MatrixX3d> &paths);

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &zeroForm,
             const Eigen::VectorXd &oneForm, const Eigen::MatrixXi &FESigns,
             const Eigen::MatrixX3d &B1, const Eigen::MatrixX3d &B2,
             const Eigen::MatrixX3d &B3,
             const double &value, const bool &updateGeometry,
             std::vector<Eigen::MatrixX3d> &paths,
             std::vector<std::vector<int>> &pathIds,
             Eigen::MatrixX3d &NV, Eigen::MatrixX3i &NF, Eigen::VectorXd &otherParameterization);
}
}
