#pragma once

#include <Eigen/Core>
#include "../utils/opengl/glfw/Viewer.h"
#include "../action/types.h"
#include "./types.h"

namespace rp {
namespace viewer {
  using MouseDownCallback = std::function<bool(rp::opengl::glfw::Viewer &viewer, int, int)>;
  using MouseMoveCallback = std::function<bool(rp::opengl::glfw::Viewer &viewer, int, int)>;
  using MouseUpCallback = std::function<bool(rp::opengl::glfw::Viewer &viewer, int, int)>;

  MouseDownCallback use_mouse_down(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                                   SelectMode &selectMode,
                                   rp::action::GeometrySetter &geometrySetter,
                                   rp::action::PathSetter &pathSetter,
                                   rp::action::PointsSetter &pointSetter,
                                   rp::action::EdgePathStateSetter &setGeodesicBrushPath,
                                   rp::action::EdgePathStateSetter &setGeodesicPrescriptionPath,
                                   rp::action::IndexExchangeStateSetter &exchangeIndices,
                                   rp::action::IndexExchangeStateSetter &annihilateIndices,
                                   rp::action::IndexExchangeStateSetter &addSingularityPair,
                                   rp::action::EdgePathBrushStateSetter &setBrushpath);

  MouseMoveCallback use_mouse_move(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                                   SelectMode &selectMode, 
                                   rp::action::EdgePathStateSetter &setGeodesicBrushPath,
                                   rp::action::EdgePathBrushStateSetter &setBrushPath,
                                   rp::action::EdgePathBrushStateSetter &setPrescriptionPath,
                                   rp::action::VectorFieldStateSetter &computeVectorfield,
                                   rp::action::ParameterizationStateSetter &parameterize);

  MouseUpCallback use_mouse_up(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                               SelectMode &selectMode, 
                               rp::action::EdgePathStateSetter &setGeodesicBrushPath,
                               rp::action::EdgePathBrushStateSetter &setBrushPath,
                               rp::action::EdgePathBrushStateSetter &setPrescriptionPath,
                               rp::action::VectorFieldStateSetter &computeVectorfield,
                               rp::action::ParameterizationStateSetter &parameterize);

  std::function<void()> use_reset_interaction();
}
}
