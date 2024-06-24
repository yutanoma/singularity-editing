#include "./parameterize.h"

namespace rp {
namespace action {
  ParameterizationStateSetter use_parameterize(AngleSetter &angleSetter, IndicesSetter &indicesSetter) {
    ParameterizationStateSetter parameterize = [&]()
    {
      rp::dispatcher::parameterize();

      Eigen::VectorXd a = rp::store::get_parameterization();
      Eigen::VectorXd da = rp::store::get_parameterization_diff();
      Eigen::MatrixXi EV = rp::store::get_ev();
      Eigen::MatrixXi FE = rp::store::get_fe();
      angleSetter(a, da, EV, FE);

      auto F = rp::store::get_f();
      auto V = rp::store::get_v();
      auto indices = rp::store::get_parameterization_indices();
      indicesSetter(F, V, indices);
    };

    return parameterize;
  };
}
}
