#pragma once

#include <Eigen/Core>
#include <functional>

namespace rp {
namespace action {
  // ジオメトリをセット
  using GeometrySetter = std::function<void(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V)>;
  // パスをセット
  using PathSetter = std::function<void(const Eigen::MatrixX3d &path, const Eigen::MatrixX3d &color)>;
  // ベクトル場をセット
  using FieldSetter = std::function<void(const Eigen::MatrixXd &vectorField, const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V)>;
  // 点をセット
  using PointsSetter = std::function<void(const Eigen::MatrixX3d &path, const Eigen::MatrixX3d &color)>;
  // 色をセット
  using ColorSetter = std::function<void(const Eigen::MatrixX3d &colors)>;
  using AngleSetter = std::function<void(const Eigen::VectorXd &data, const Eigen::VectorXd &ddata, const Eigen::MatrixXi &EV, const Eigen::MatrixXi &FE)>;
  // parameterizationをセット
  using ParameterSetter = std::function<void(const Eigen::VectorXd &parameterization)>;
  using IndicesSetter = std::function<void(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V, const Eigen::VectorXd &)>;

  using VectorFieldStateSetter = std::function<void()>;
  using IndexExchangeStateSetter = std::function<void(const int &fid1, const int &fid2)>;
  using PathStateSetter = std::function<void(const int &nextVid, const bool &isNewPath, Eigen::MatrixX3d &newPaths)>;
  using ParameterizationStateSetter = std::function<void()>;
  using RescalingStateSetter = std::function<void()>;
  using StripePatternsStateSetter = std::function<void()>;
  using InitializeStateSetter = std::function<void(Eigen::MatrixX3i &F, Eigen::MatrixX3d &V)>;
  using EdgePathStateSetter = std::function<void(const int &prevVid, const int &nextVid)>;
  using EdgePathBrushStateSetter = std::function<void(const int &vid0, const int &vid1, const double &ratio)>;
  using Saver = std::function<void()>;
};
};
