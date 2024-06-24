#include "./open_path.h"

#include <igl/remove_unreferenced.h>

#include "./remove_duplicate_vertices.h"
#include "./divide_borders.h"

namespace rp {
void open_path(const std::vector<Eigen::MatrixX3d> &paths, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, std::vector<std::vector<int>> &pathVids) {
  std::vector<std::vector<std::array<int, 2>>> vidsPerFace;
  rp::cut_mesh_with_path::process(F, V, paths, vidsPerFace, pathVids);
  rp::subdivide_by_path::process(V, F, vidsPerFace);
  rp::divide_borders::process(V, F, pathVids);

  Eigen::MatrixX3i NF;
  Eigen::MatrixX3d NV;
  Eigen::MatrixXi I;
  igl::remove_unreferenced(V, F, NV, NF, I);

  F = NF;
  V = NV;
}
}
