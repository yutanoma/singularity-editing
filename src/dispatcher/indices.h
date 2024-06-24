#pragma once

#include <Eigen/Core>
#include "../store/geometry.h"

namespace rp {
namespace dispatcher {
  void prescribe_parameterization_indices(const Eigen::VectorXd &indices);
}
}
