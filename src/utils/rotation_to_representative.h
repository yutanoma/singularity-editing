// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <igl/edge_topology.h>
#include <igl/gaussian_curvature.h>
#include <igl/igl_inline.h>
#include <igl/local_basis.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <vector>

namespace rp {
// Converts a rotation-angle representation + global rotation on face 0, to a
// representative format encoding an N-RoSy. Note: the rotation angles must
// respect the triviality conditions within the degree N, or the results are
// only least-square fitting. Inputs:
//  V:              #V x 3 vertex coordinates
//  F:              #F x 3 face vertex indices
//  EV:             #E x 2 edges to vertices indices
//  EF:             #E X 2 edges to faces indices
//  B1, B2:         #F x 3 matrices representing the local base of each face.
//  rotationAngles: #E angles that encode the rotation angles (deviation from
//  parallel transport) N:              the degree of the field globalRotation:
//  The angle between the vector on the first face and its basis in radians
// Outputs:
//  representative: #F x 3 representative vectors on the faces.

void rotation_to_representative(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
    const Eigen::MatrixXd &B1, const Eigen::MatrixXd &B2,
    const Eigen::VectorXd &rotationAngles, const int N,
    const double globalRotation, Eigen::MatrixXd &representative);

// Version without local basis as input
void rotation_to_representative(const Eigen::MatrixXd &V,
                                const Eigen::MatrixXi &F,
                                const Eigen::MatrixXi &EV,
                                const Eigen::MatrixXi &EF,
                                const Eigen::VectorXd &rotationAngles,
                                const int N, double globalRotation,
                                Eigen::MatrixXd &representative);
}  // namespace rp
