#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include "../action/types.h"

namespace rp {
namespace viewer {
rp::action::FieldSetter use_vector_field_setter(rp::opengl::glfw::Viewer &viewer);

std::function<void(void)> use_vector_field_clearer(rp::opengl::glfw::Viewer &viewer);
}
}
