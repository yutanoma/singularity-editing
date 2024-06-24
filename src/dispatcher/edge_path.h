#pragma once

#include <Eigen/Core>
#include "../store/geometry.h"

namespace rp {
namespace dispatcher {
  void set_brush_paths(const int &prevVid, const int &nextVid, Eigen::MatrixX3d &newPaths);

  void set_brush_stroke(const int &_v0, const int &_v1, const double &_ratio, Eigen::MatrixX3d &newPaths);

  void set_index_prescription_paths();

  void set_prescription_path_index();
};
};
