#pragma once

#include <Eigen/Core>
#include "vector_utils.h"

namespace rp {
  void rotate_vector_field(const Eigen::VectorXd &originalVectorAngles, const double &angle,
                           Eigen::VectorXd &rotatedVectorAngles);
}
