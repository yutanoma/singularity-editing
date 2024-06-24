// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "index_prescription.h"

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
    Eigen::VectorXd &rotationAngles, double &linfError)
{
  index_prescription(V, F, EV, EF, innerEdges, basisCycles, cycleCurvature, cycleIndices, ldltSolver, N, weightMatrix, true, rotationAngles, linfError);
}

void index_prescription(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, const Eigen::VectorXi &innerEdges,
    const Eigen::SparseMatrix<double> &basisCycles,
    const Eigen::VectorXd &cycleCurvature, const Eigen::VectorXd &cycleIndices,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &ldltSolver,
    const int N, const Eigen::SparseMatrix<double> &weightMatrix, 
    const bool needsRefactorize, Eigen::VectorXd &rotationAngles, double &linfError)
{
  using namespace Eigen;
  using namespace std;

  VectorXd cycleNewCurvature =
      cycleIndices.cast<double>() * (2.0 * igl::PI / (double)N);

  SparseMatrix<double> AD = basisCycles * weightMatrix;
  
  // for (int i = 0; i < basisCycles.rows(); i++) {
  //   bool ok = false;
  //   for (int j = 0; j < basisCycles.cols(); j++) {
  //     if (AD.coeff(i, j) != 0) {
  //       ok = true;
  //       break;
  //     }
  //   }
  //   if (!ok) {
  //     std::cout << "row: " << i << std::endl;
  //   }
  // }

  // for (int i = 0; i < basisCycles.cols(); i++) {
  //   bool ok = false;
  //   for (int j = 0; j < basisCycles.rows(); j++) {
  //     if (AD.coeff(j, i) != 0) {
  //       ok = true;
  //       break;
  //     }
  //   }
  //   if (!ok) {
  //     std::cout << "col: " << i << std::endl;
  //   }
  // }

  // Initialize solver if never before
  if (!ldltSolver.rows() || needsRefactorize) {
    SparseMatrix<double> AAt = AD * AD.transpose();
    ldltSolver.compute(AAt);
  }

  std::cout << "solver computed: " << ldltSolver.info() << std::endl;

  std::cout << "bsc: " << basisCycles.rows()
            << ", cc: " << cycleCurvature.rows()
            << ", cnc: " << cycleNewCurvature.rows()
            << ", weightMatrix: " << weightMatrix.rows()
            << ", innerEdges: " << innerEdges.rows() << std::endl;

  assert(ldltSolver.info() == Eigen::Success);

  Eigen::VectorXd solvedVec = ldltSolver.solve(-cycleCurvature + cycleNewCurvature);

  auto mat = weightMatrix * weightMatrix.transpose() * basisCycles.transpose();

  VectorXd innerRotationAngles = weightMatrix * weightMatrix.transpose() * basisCycles.transpose() * solvedVec;

  rotationAngles.resize(EV.rows());
  rotationAngles.setZero();

  for (int i = 0; i < innerEdges.rows(); i++)
    rotationAngles(innerEdges(i)) = innerRotationAngles(i);

  auto resVect = basisCycles * innerRotationAngles;
  auto diffVect = (resVect - (-cycleCurvature + cycleNewCurvature));
  linfError = diffVect.lpNorm<Infinity>();

  // if (linfError > 1) {
  //   std::cout << diffVect << std::endl;
  // }

  std::cout << "index_prescription linfError: " << linfError << ", rotationAngles norm" << rotationAngles.norm() << std::endl;
  // for (int i = diffVect.rows() - 100; i < diffVect.rows(); i++) {
  //   std::cout << diffVect(i) << ", " << cycleCurvature(i) << ", " << cycleNewCurvature(i) << std::endl;
  // }
}
}  // namespace rp
