#pragma once

#include <igl/boundary_loop.h>
#include <igl/remove_unreferenced.h>
#include <iostream>

#include <Eigen/Core>
#include <array>
#include <vector>

namespace rp {
namespace divide_borders {
  void process(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
               // 切り開く点のindices。ループの場合は最後の要素と最初の要素が一致するようにする。
               const std::vector<std::vector<int>> &positions);

  void process(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
               // 切り開く点のindices
               const std::vector<Eigen::VectorXi> &positions,
               std::vector<std::vector<int>> &v_connection);

  void process(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
               // 切り開く点のindices
               const std::vector<Eigen::VectorXi> &positions,
               std::vector<std::vector<int>> &v_connection,
               std::vector<std::array<int, 2>> &new_vertices);
};  // namespace divide_borders
};  // namespace rp
