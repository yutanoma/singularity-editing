#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include "../action/types.h"

namespace rp {
namespace viewer {
rp::action::GeometrySetter use_geometry_setter(rp::opengl::glfw::Viewer &viewer, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, int viewerId);
}
}
