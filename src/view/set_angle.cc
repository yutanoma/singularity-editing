#include "./set_angle.h"
#include <igl/writeOBJ.h>
#include <igl/isolines_map.h>
#include <igl/opengl/destroy_shader_program.h>
#include <igl/opengl/create_shader_program.h>
#include "../../data/path.h"
#include "../utils/vector_utils.h"
#include "../utils/angle_color.h"

namespace rp {
namespace viewer {
rp::action::AngleSetter use_angle_setter(rp::opengl::glfw::Viewer &viewer) {
  auto parameterSetter = [&](const Eigen::VectorXd &data, const Eigen::VectorXd &ddata,
                             const Eigen::MatrixXi &EV, const Eigen::MatrixXi &FE)
  {
    viewer.data_list[0].set_angles(data, ddata, EV, FE);
  };

  return parameterSetter;
};
}
}
