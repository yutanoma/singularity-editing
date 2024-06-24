#include "./button_operations.h"

namespace rp {
namespace dispatcher {
  void compute_rescalings() {
    rp::store::compute_rescalings();
  }

  void compute_stripe_patterns() {
    rp::store::set_default_stripe_patterns();
  };

  void parameterize() {
    rp::store::recompute_parameterization();
  }

  void reset_parameterization() {
    rp::store::temporary_paths::reset_paths();
    rp::store::reset_parameterization();
  }

  void reset_rescalings() {
    rp::store::reset_rescalings();
  }
};
};
