// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>
#include <igl/gaussian_curvature.h>
#include <igl/igl_inline.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <iostream>

#include <Eigen/Core>
#include <cmath>
#include <vector>

#include "vector_utils.h"

namespace rp {
// Computes the dual-edge-based rotation angles that are required to reproduce a
// prescribed set of indices on the dual cycles of the mesh. In case the sum of
// curvature is not consistent with the topology, the system is solved in least
// squares and unexpected singularities may appear elsewhere. linfError will
// mostl like be far from zero. Inputs:
//  V:          #V by 3 vertex coordinates
//  F:          #F by 3 face vertex indices
//  EV:         #E by 3 edges
//  innerEdges: #iE the subset from EV of inner (non-boundary) edges.
//  basisCycles:#c X #E the basis cycles matrix (obtained from
//  directional::dual_cycles indices:    #c the prescribed index around each
//  cycle. They should add up to N*Euler_characteristic of the mesh.
//  cycleCurvature: #c the original curvature for each basis cycle.
//  solver: The Simplicial LDLT solver used to solver the problem. It will only
//  prefactor the matrix once upon the first call to the function; the state of
//  the solver solely depends on the basisCycles, therefore it only needs to be
//  reset if the basisCycles matrix changed. N: the degree of the field.
// Output:
//  rotationAngles: #iE rotation angles (difference from parallel transport) per
//  inner dual edge linfError: l_infinity error of the computation. If this is
//  not approximately 0, the prescribed indices are likely inconsistent (don't
//  add up to the correct sum).
void index_prescription(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, const Eigen::VectorXi &innerEdges,
    const Eigen::SparseMatrix<double> &basisCycles,
    const Eigen::VectorXd &cycleCurvature, const Eigen::VectorXd &cycleIndices,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &ldltSolver,
    const int N, const Eigen::SparseMatrix<double> &weightMatrix,
    Eigen::VectorXd &rotationAngles, double &linfError);

void index_prescription(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, const Eigen::VectorXi &innerEdges,
    const Eigen::SparseMatrix<double> &basisCycles,
    const Eigen::VectorXd &cycleCurvature, const Eigen::VectorXd &cycleIndices,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &ldltSolver,
    const int N, const Eigen::SparseMatrix<double> &weightMatrix, 
    const bool needsRefactorize, Eigen::VectorXd &rotationAngles, double &linfError);
}  // namespace rp
