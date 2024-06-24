#include "./set_path.h"
#include "../utils/line_cylinders.h"
#include <directional/point_spheres.h>

namespace rp {
namespace viewer {
namespace {
  int sphereMeshId = -1;
  Eigen::MatrixX3d spherePoints(0, 3);
  Eigen::MatrixX3d sphereColors(0, 3);

  int pathsMeshId = -1;
  Eigen::MatrixX3d P1(0, 3), P2(0, 3);
  Eigen::MatrixX3d pathsColors(0, 3);
};

rp::action::PointsSetter use_point_setter(rp::opengl::glfw::Viewer &viewer) {
  auto pointSetter = [&](const Eigen::MatrixX3d &points, const Eigen::MatrixX3d &color)
  {
    if (sphereMeshId == -1) {
      sphereMeshId = viewer.append_mesh();
    }

    int initialSphereRows = spherePoints.rows();
    spherePoints.conservativeResize(initialSphereRows + points.rows(), 3);
    sphereColors.conservativeResize(initialSphereRows + points.rows(), 3);

    for (int i = 0; i < points.rows(); i++) {
      spherePoints.row(i + initialSphereRows) = points.row(i);

      if (color.rows() > i) {
        sphereColors.row(i + initialSphereRows) = color.row(i);
      } else if (color.rows() > 0) {
        sphereColors.row(i + initialSphereRows) = color.row(0);
      } else {
        sphereColors.row(i + initialSphereRows) << 233. / 255., 76. / 255., 38. / 255.;
      }
    }

    Eigen::MatrixXd sV, C;
    Eigen::MatrixXi sF;

    directional::point_spheres(spherePoints, 0.5, sphereColors, 12, sV, sF, C);

    viewer.data_list[sphereMeshId].clear();
    viewer.data_list[sphereMeshId].set_mesh(sV, sF);
    viewer.data_list[sphereMeshId].set_colors(C);
    viewer.data_list[sphereMeshId].show_lines = false;
  };

  return pointSetter;
}

rp::action::PathSetter use_path_setter(rp::opengl::glfw::Viewer &viewer) {
  auto pathSetter = [&](const Eigen::MatrixX3d &path, const Eigen::MatrixX3d &color)
  {
    if (path.rows() == 0) {
      return;
    }

    if (pathsMeshId == -1) {
      pathsMeshId = viewer.append_mesh();
    }

    int initialPathsSize = P1.rows();

    P1.conservativeResize(initialPathsSize + path.rows() - 1, 3);
    P2.conservativeResize(initialPathsSize + path.rows() - 1, 3);
    pathsColors.conservativeResize(initialPathsSize + path.rows() - 1, 3);

    for (int i = 0; i < path.rows() - 1; i++) {
      P1.row(initialPathsSize + i) = path.row(i);
      P2.row(initialPathsSize + i) = path.row(i + 1);

      if (color.rows() > i) {
        pathsColors.row(i + initialPathsSize) = color.row(i);
      } else if (color.rows() > 0) {
        pathsColors.row(i + initialPathsSize) = color.row(0);
      } else {
        pathsColors.row(i + initialPathsSize) << 233. / 255., 76. / 255., 38. / 255.;
      }
    }

    Eigen::MatrixXd V, C;
    Eigen::MatrixXi F;

    rp::line_cylinders(P1, P2, 0.3, pathsColors, 4, V, F, C);

    viewer.data_list[pathsMeshId].clear();
    viewer.data_list[pathsMeshId].set_mesh(V, F);
    viewer.data_list[pathsMeshId].set_colors(C);
    viewer.data_list[pathsMeshId].show_lines = false;
  };

  return pathSetter;
};

std::function<void(void)> use_path_clearer(rp::opengl::glfw::Viewer &viewer) {
  auto clearer = [&]()
  {
    viewer.data_list[sphereMeshId].clear();
    viewer.data_list[pathsMeshId].clear();
  };

  return clearer;
};
};
};
