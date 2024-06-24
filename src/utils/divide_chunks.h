#pragma once

#include <igl/boundary_loop.h>
#include <igl/remove_unreferenced.h>
#include <iostream>

#include <Eigen/Core>
#include <vector>

// 全ての面に対してchunkに分割する
// 接続状況から判定する

namespace rp {
namespace divide_chunks {

void process(
    Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
    std::vector<std::vector<int>> &v_connection,
    std::vector<Eigen::MatrixX3i> &result_faces,
    std::vector<Eigen::MatrixX3d> &result_vertices,
    std::vector<std::vector<std::vector<int>>> &result_border,
    std::vector<std::vector<std::vector<std::vector<int>>>> &result_connection,
    std::vector<std::vector<int>> &initial_boundary_loops);
};
};  // namespace rp
