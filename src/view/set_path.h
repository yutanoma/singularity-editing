#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include "../action/types.h"

namespace rp {
namespace viewer {
rp::action::PointsSetter use_point_setter(rp::opengl::glfw::Viewer &viewer);

rp::action::PathSetter use_path_setter(rp::opengl::glfw::Viewer &viewer);

std::function<void(void)> use_path_clearer(rp::opengl::glfw::Viewer &viewer);
}
}
