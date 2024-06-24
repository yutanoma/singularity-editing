#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

namespace rp {
namespace cut_mesh_with_path {
void process(const Eigen::MatrixX3i &F, Eigen::MatrixX3d &V, const std::vector<Eigen::MatrixX3d> &paths, std::vector<std::vector<std::array<int, 2>>> &vidsPerFace, std::vector<std::vector<int>> &pathVids);
}
}
