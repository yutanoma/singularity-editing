#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include "../action/types.h"

namespace rp {
namespace viewer {
  rp::action::AngleSetter use_angle_setter(rp::opengl::glfw::Viewer &viewer);
}
}
