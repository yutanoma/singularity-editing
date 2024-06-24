#include "set_vector_field.h"

#include <directional/glyph_lines_raw.h>
#include <directional/visualization_schemes.h>

namespace rp {
namespace viewer {
namespace {
  int meshId = -1;
}

rp::action::FieldSetter use_vector_field_setter(rp::opengl::glfw::Viewer &viewer) {
  auto vectorFieldSetter = [&](const Eigen::MatrixXd &vectorField, const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V)
  {
    Eigen::MatrixXi FField, _F = F;
    Eigen::MatrixXd VField, CField, _V = V;

    if (meshId == -1) {
      meshId = viewer.append_mesh();
    }

    Eigen::MatrixXd colors(1, 3);
    colors << 1, 1, 1;

    directional::glyph_lines_raw(_V, _F, vectorField, colors, VField, FField, CField, 5);

      std::cout << "meshId: " << meshId << std::endl;

    std::cout << VField.rows() << ", " << FField.rows() << ", " << V.rows() << ", " << F.rows() << std::endl;

    viewer.data_list[meshId].set_mesh(VField, FField);
    viewer.data_list[meshId].set_colors(CField);
    viewer.data_list[meshId].show_lines = false;
    viewer.data_list[meshId].show_faces = true;
  };

  return vectorFieldSetter;
};

std::function<void(void)> use_vector_field_clearer(rp::opengl::glfw::Viewer &viewer) {
  auto vectorFieldClearer = [&]()
  {
    if (meshId != -1)
    {
      viewer.data_list[meshId].clear();
      meshId = -1;
    }
  };

  return vectorFieldClearer;
};
}
}
