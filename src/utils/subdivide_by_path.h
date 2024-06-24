#pragma once

#include <eigen/Core>
#include <vector>
#include <array>
#include <iostream>

namespace rp {
namespace subdivide_by_path {
// 各faceの上に乗っているパスのvidを指定すると、そのパスでメッシュをsubdivideする
void process(const Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, std::vector<std::vector<std::array<int, 2>>> &pathOnFaces);
}
}
