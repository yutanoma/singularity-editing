#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include "../action/types.h"

namespace rp {
namespace viewer {
  rp::action::IndicesSetter use_indices_setter(rp::opengl::glfw::Viewer &viewer, double radius);
}
}
