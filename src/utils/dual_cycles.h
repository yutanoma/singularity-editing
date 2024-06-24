// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <directional/tree.h>
#include <igl/boundary_loop.h>
#include <igl/colon.h>
#include <igl/edge_topology.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <iostream>

#include <Eigen/Core>
#include <unordered_map>
#include <vector>
#include <stack>

namespace rp {
// Creates the set of independent dual cycles (closed loops of connected faces
// that cannot be morphed to each other) on a mesh. Primarily used for index
// prescription. The basis cycle matrix first contains #V-#b cycles for every
// inner vertex (by order), then #b boundary cycles, and finally 2*g generator
// cycles around all handles. Total #c cycles.The cycle matrix sums information
// on the dual edges between the faces, and is indexed into the inner edges
// alone (excluding boundary)
// input:
//  V: #V by 3 vertices.
//  F: #F by 3 triangles.
//  EV: #E by 2 matrix of edges (vertex indices)
//  EF: #E by 2 matrix of oriented adjacent faces
// hasConstraint: #F vector if the face has constraints. We note that to enable
// constraints, hasConstraint[0] must be true. constraint: #F vector constrained
// angles of each faces. the constraint is only valid when the corresponding
// value of hasConstraint is true B1: #F by 3 the basic direction of each face
// output:
//  basisCycles:    #c by #iE basis cycles
//  cycleCurvature:   #c by 1 curvatures of each cycle (for inner-vertex cycles,
//  simply the Gaussian curvature. vertex2cycle:     #v by 1 map between vertex
//  and corresponding cycle (for comfort of input from the user's side; inner
//  vertices map to their cycles, boundary vertices to the bigger boundary
//  cycle. innerEdges:       #iE by 1 the subset of #EV that are inner edges,
//  and with the same ordering as the columns of basisCycles.

void dual_cycles(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                 const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF,
                 Eigen::MatrixXi& FE, const std::vector<bool>& hasConstraint,
                 const std::vector<double>& constraint,
                 const Eigen::MatrixX3d& B1,
                 Eigen::SparseMatrix<double>& basisCycles,
                 Eigen::VectorXd& cycleCurvature, Eigen::VectorXi& vertex2cycle,
                 Eigen::VectorXi& innerEdges, int& innerVerticesSize);
}  // namespace rp
