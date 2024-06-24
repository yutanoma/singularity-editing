#include "./set_geometry.h"
#include <igl/writeOBJ.h>
#include "../../data/path.h"

namespace rp {
namespace viewer {
rp::action::GeometrySetter use_geometry_setter(rp::opengl::glfw::Viewer &viewer, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, int viewerId) {
  auto geometrySetter = [&](const Eigen::MatrixX3i &_F, const Eigen::MatrixX3d &_V)
  {
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(_V, _F);
    viewer.data_list[0].show_lines = false;

    Eigen::MatrixX3d colors(_F.rows(), 3);
    for (int i = 0; i < colors.rows(); i++) {
      colors.row(i) << 235. / 255., 242. / 255., 244. / 255.;
    }

    viewer.data_list[0].set_colors(colors);
    V = _V;
    F = _F;
    // igl::writeOBJ(DATA_PATH "/models/bunny_opened.obj", V, F);
  };

  return geometrySetter;
};
}
}
