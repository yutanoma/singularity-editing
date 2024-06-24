#include "./set_geometry.h"
#include <igl/writeOBJ.h>
#include "../../data/path.h"
#include "../utils/vector_utils.h"
#include "../utils/angle_color.h"

namespace rp {
namespace viewer {
rp::action::ColorSetter use_color_setter(rp::opengl::glfw::Viewer &viewer) {
  auto colorSetter = [&](const Eigen::MatrixX3d &colors)
  {
    viewer.data_list[0].set_colors(colors);
  };

  return colorSetter;
};

rp::action::ParameterSetter use_parameter_setter(rp::opengl::glfw::Viewer &viewer) {
  auto parameterSetter = [&](const Eigen::VectorXd &data)
  {
    viewer.data_list[0].set_data(data);
  };

  return parameterSetter;
};
}
}
