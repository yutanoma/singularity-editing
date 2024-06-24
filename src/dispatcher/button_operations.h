#pragma once

#include <Eigen/Core>
#include "../store/geometry.h"
#include "../store/temporary_paths.h"

namespace rp {
namespace dispatcher {
  void compute_rescalings();

  void compute_stripe_patterns();

  void parameterize();

  void reset_parameterization();

  void reset_rescalings();
};
};
