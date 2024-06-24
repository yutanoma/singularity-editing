#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include "../action/types.h"

namespace rp {
namespace viewer {
  rp::action::ColorSetter use_color_setter(rp::opengl::glfw::Viewer &viewer);

  rp::action::ParameterSetter use_parameter_setter(rp::opengl::glfw::Viewer &viewer);
}
}
