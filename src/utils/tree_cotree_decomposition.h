#pragma once

#include <Eigen/Core>
#include <directional/tree.h>
#include <igl/boundary_loop.h>
#include <igl/local_basis.h>
#include <igl/gaussian_curvature.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <igl/edge_topology.h>

namespace rp {
  void tree_cotree_decomposition(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, const int &rowOffset, std::vector<Eigen::Triplet<double>> &basisCycleTriplets);
}
