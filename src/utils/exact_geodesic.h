// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <igl/PI.h>
#include <iostream>

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <vector>
#include <array>

namespace rp {
namespace exact_geodesic {
// Exact geodesic algorithm for triangular mesh with the implementation from
// https://code.google.com/archive/p/geodesic/, and the algorithm first
// described by Mitchell, Mount and Papadimitriou in 1987
//
// Inputs:
//   V  #V by 3 list of 3D vertex positions
//   F  #F by 3 list of mesh faces
//   VS #VS by 1 vector specifying indices of source vertices
//   FS #FS by 1 vector specifying indices of source faces
//   VT #VT by 1 vector specifying indices of target vertices
//   FT #FT by 1 vector specifying indices of target faces
// Output:
//   D  #VT+#FT by 1 vector of geodesic distances of each target w.r.t. the
//   nearest one in the source set
//
// Note:
//      Specifying a face as target/source means its center.
//
template <typename DerivedV, typename DerivedF, typename DerivedVS,
          typename DerivedFS, typename DerivedVT, typename DerivedFT,
          typename DerivedD>
inline void exact_geodesic(const Eigen::MatrixBase<DerivedV> &V,
                           const Eigen::MatrixBase<DerivedF> &F,
                           const Eigen::MatrixBase<DerivedVS> &VS,
                           const Eigen::MatrixBase<DerivedFS> &FS,
                           const Eigen::MatrixBase<DerivedVT> &VT,
                           const Eigen::MatrixBase<DerivedFT> &FT,
                           Eigen::PlainObjectBase<DerivedD> &D);

// 点p1から点p2までの間をedgeに沿って通った時の最小の長さ（と思われる）のパス
void get_crack_path(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, const int p1,
                    const int p2, std::vector<int> &result_path);

void split_path(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                // 切り開く点のindices
                const std::vector<Eigen::VectorXi> &positions,
                // 切り開いた箇所の座標
                std::vector<Eigen::MatrixX3d> &result_points);

void split_path(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                // 切り開く点のindices
                const std::vector<Eigen::VectorXi> &positions,
                // verticesやfacesを切り開くかどうか
                const bool update_geometry,
                std::vector<Eigen::VectorXi> &result_indices);

// 点の位置からパスを作る。upadte_geometryをtrueにすると当該のパスで三角形を分割する。
void split_path(
    Eigen::MatrixX3d& V, Eigen::MatrixX3i& F,
    // 切り開く点のindices。途中で交差するものは禁止、必ず既存の点から分岐するようにする
    const std::vector<Eigen::VectorXi>& positions,
    // verticesやfacesを切り開くかどうか
    const bool update_geometry,
    // 切り開いた箇所の座標
    std::vector<Eigen::MatrixX3d>& result_points,
    // 切り開いた箇所のindices、update_geometryがtrueだった時にのみ編集する
    std::vector<Eigen::VectorXi>& result_indices,
    std::vector<std::vector<std::pair<std::array<int, 2>, double>>> &edgePaths);
}  // namespace exact_geodesic
}  // namespace rp
