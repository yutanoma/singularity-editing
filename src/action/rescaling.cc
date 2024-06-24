#include "./rescaling.h"

namespace rp {
namespace action {
RescalingStateSetter use_rescaling_setter(AngleSetter &angleSetter) {
  RescalingStateSetter rescalingSetter = [&]()
  {
    rp::dispatcher::compute_rescalings();

    Eigen::VectorXd a = rp::store::get_parameterization();
    Eigen::VectorXd da = rp::store::get_parameterization_diff();
    Eigen::MatrixXi EV = rp::store::get_ev();
    Eigen::MatrixXi FE = rp::store::get_fe();
    angleSetter(a, da, EV, FE);
  };

  return rescalingSetter;
}

RescalingStateSetter use_rescaling_resetter(AngleSetter &angleSetter) {
  RescalingStateSetter rescalingResetter = [&]()
  {
    rp::dispatcher::reset_rescalings();

    Eigen::VectorXd a = rp::store::get_parameterization();
    Eigen::VectorXd da = rp::store::get_parameterization_diff();
    Eigen::MatrixXi EV = rp::store::get_ev();
    Eigen::MatrixXi FE = rp::store::get_fe();
    angleSetter(a, da, EV, FE);
  };

  return rescalingResetter;
}
}
}

