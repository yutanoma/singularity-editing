#pragma once

#include <Eigen/Core>
#include "../utils/opengl/glfw/Viewer.h"
#include "./types.h"
#include "../dispatcher/button_operations.h"
#include "../store/geometry.h"

namespace rp {
namespace action {
  ParameterizationStateSetter use_parameterize(AngleSetter &angleSetter, IndicesSetter &indicesSetter);
}
}
