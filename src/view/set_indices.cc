#include "./set_indices.h"

#include "../utils/singularity_spheres.h"

namespace rp {
namespace viewer {
namespace {
  int meshId = -1;
  double _radius = 0;
}

rp::action::IndicesSetter use_indices_setter(rp::opengl::glfw::Viewer &viewer, double radius) {
  _radius = radius;

  auto indicesSetter = [&](const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V, const Eigen::VectorXd &indices) {
    if (meshId == -1) {
      meshId = viewer.append_mesh();
    }

    Eigen::MatrixXd sV, sC;
    Eigen::MatrixXi sF;
    rp::singularity_spheres(F, V, indices, sF, sV, sC, _radius);

    viewer.data_list[meshId].clear();
    viewer.data_list[meshId].set_mesh(sV, sF);
    viewer.data_list[meshId].set_colors(sC);
    viewer.data_list[meshId].show_lines = false;
  };

  return indicesSetter;
}
}
}
