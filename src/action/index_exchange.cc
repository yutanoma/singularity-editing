#include "./index_exchange.h"

namespace rp {
namespace action {
// exchange the indices of the selected two faces
IndexExchangeStateSetter use_index_exchange_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter) {
  IndexExchangeStateSetter indexExchangeSetter = [&](const int &fid1, const int &fid2)
  {
    auto indices = rp::store::get_parameterization_indices();
    std::swap(indices(fid1), indices(fid2));

    rp::dispatcher::prescribe_parameterization_indices(indices);

    Eigen::VectorXd a = rp::store::get_parameterization();
    Eigen::VectorXd da = rp::store::get_parameterization_diff();
    Eigen::MatrixXi EV = rp::store::get_ev();
    Eigen::MatrixXi FE = rp::store::get_fe();
    angleSetter(a, da, EV, FE);

    auto F = rp::store::get_f();
    auto V = rp::store::get_v();
    indicesSetter(F, V, indices);
  };

  return indexExchangeSetter;
};

// set the indices of the selected two faces to +1 and -1
IndexExchangeStateSetter use_add_singularity_pair_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter) {
  IndexExchangeStateSetter singularityPairSetter = [&](const int &fid1, const int &fid2) {
    auto indices = rp::store::get_parameterization_indices();

    if (indices(fid1) == 0 && indices(fid2) == 0) {
      indices(fid1) = +1;
      indices(fid2) = -1;

      rp::dispatcher::prescribe_parameterization_indices(indices);

      Eigen::VectorXd a = rp::store::get_parameterization();
      Eigen::VectorXd da = rp::store::get_parameterization_diff();
      Eigen::MatrixXi EV = rp::store::get_ev();
      Eigen::MatrixXi FE = rp::store::get_fe();
      angleSetter(a, da, EV, FE);

      auto F = rp::store::get_f();
      auto V = rp::store::get_v();
      indicesSetter(F, V, indices);
    }
  };

  return singularityPairSetter;
}

// set the indices of the selected two faces to 0
IndexExchangeStateSetter use_index_annihilation_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter) {
  IndexExchangeStateSetter indexExchangeSetter = [&](const int &fid1, const int &fid2)
  {
    auto indices = rp::store::get_parameterization_indices();

    if (indices(fid1) + indices(fid2) == 0) {
      indices(fid1) = 0;
      indices(fid2) = 0;

      rp::dispatcher::prescribe_parameterization_indices(indices);

      Eigen::VectorXd a = rp::store::get_parameterization();
      Eigen::VectorXd da = rp::store::get_parameterization_diff();
      Eigen::MatrixXi EV = rp::store::get_ev();
      Eigen::MatrixXi FE = rp::store::get_fe();
      angleSetter(a, da, EV, FE);

      auto F = rp::store::get_f();
      auto V = rp::store::get_v();
      indicesSetter(F, V, indices);
    }
  };

  return indexExchangeSetter;
};
}
}
