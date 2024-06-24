#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

namespace rp {
  bool line_cylinders(const Eigen::MatrixXd &P1,
                      const Eigen::MatrixXd &P2,
                      const double &radius,
                      const Eigen::MatrixXd &cyndColors,
                      const int res,
                      Eigen::MatrixXd &V,
                      Eigen::MatrixXi &T,
                      Eigen::MatrixXd &C);
}
