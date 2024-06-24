#include "./indices.h"

namespace rp {
namespace dispatcher {
  void prescribe_parameterization_indices(const Eigen::VectorXd &indices){
    rp::store::prescribe_parameterization_indices(indices);
  };
}
}
