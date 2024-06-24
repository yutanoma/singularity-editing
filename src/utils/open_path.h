#pragma once

#include <Eigen/Core>
#include <vector>
#include <array>

#include "cut_mesh_with_path.h"
#include "subdivide_by_path.h"

namespace rp {
  void open_path(const std::vector<Eigen::MatrixX3d> &paths, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, std::vector<std::vector<int>> &pathVids);
}