#include "soft_constrained_index_prescription.h"

#include <iostream>

namespace rp {
  void soft_constrained_index_prescription(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                           const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                                           const Eigen::VectorXi &innerEdges,
                                           const Eigen::VectorXd &angleJumpsPerEdge,
                                           const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                                           const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                                           const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                                           std::vector<Eigen::Triplet<double>> &hardConstraintMatrixTriplets,
                                           const std::vector<double> &hardConstraintVals,
                                           std::vector<Eigen::Triplet<double>> &softConstraintMatrixTriplets,
                                           const std::vector<double> &softConstraintVals,
                                           Eigen::SparseMatrix<double> &matrix,
                                           Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver,
                                           Eigen::VectorXd &cycleCurvature,
                                           Eigen::VectorXd &cycleIndices) {
    // 1. まずbasisCyclesMatrixとhardConstraintを加えたものについてfactorizeしてindex_prescriptionする
    for (int i = 0; i < hardConstraintMatrixTriplets.size(); i++) {
      int row = hardConstraintMatrixTriplets[i].row();
      int col = hardConstraintMatrixTriplets[i].col();
      double val = hardConstraintMatrixTriplets[i].value();
      hardConstraintMatrixTriplets[i] = Eigen::Triplet<double>(row + basisCyclesMatrix.rows(), col, val);
      // std::cout << row << ", " << col << ", " << val << std::endl;
    }

    hardConstraintMatrixTriplets.insert(hardConstraintMatrixTriplets.begin(), basisLoopsTriplets.begin(), basisLoopsTriplets.end());

    Eigen::SparseMatrix<double> withHardConstraintMatrix(basisCyclesMatrix.rows() + hardConstraintVals.size(), basisCyclesMatrix.cols());
    withHardConstraintMatrix.setFromTriplets(hardConstraintMatrixTriplets.begin(), hardConstraintMatrixTriplets.end());

    Eigen::VectorXd actualCurvatureHardConstraints = withHardConstraintMatrix * angleJumpsPerEdge;
    cycleCurvature.conservativeResize(basisCyclesMatrix.rows() + hardConstraintVals.size());
    cycleIndices.conservativeResize(basisCyclesMatrix.rows() + hardConstraintVals.size());
    for (int i = 0; i < hardConstraintVals.size(); i++) {
      int idx = basisCyclesMatrix.rows() + i;
      cycleCurvature(idx) = rp::round_pi(- hardConstraintVals[i] + actualCurvatureHardConstraints(idx));
      cycleIndices(idx) = 0;
    }

    std::cout << cycleCurvature.minCoeff() << ", " << cycleCurvature.maxCoeff() << ", " << cycleCurvature.rows() << std::endl;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> withHardConstraintsSolver;
    Eigen::VectorXd edgeJumps;
    double linfError;

    rp::index_prescription(V, F, EV, EF, innerEdges, withHardConstraintMatrix, cycleCurvature, cycleIndices, withHardConstraintsSolver, N, weightMatrix, edgeJumps, linfError);

    // 2. 次にsoftConstraintについて一番近い値に丸めてfactorizeする
    for (int i = 0; i < softConstraintMatrixTriplets.size(); i++) {
      int row = softConstraintMatrixTriplets[i].row();
      int col = softConstraintMatrixTriplets[i].col();
      double val = softConstraintMatrixTriplets[i].value();
      softConstraintMatrixTriplets[i] = Eigen::Triplet<double>(row + basisCyclesMatrix.rows() + hardConstraintVals.size(), col, val);
    }

    softConstraintMatrixTriplets.insert(softConstraintMatrixTriplets.begin(), hardConstraintMatrixTriplets.begin(), hardConstraintMatrixTriplets.end());

    Eigen::SparseMatrix<double> withAllConstraintMatrix(basisCyclesMatrix.rows() + hardConstraintVals.size() + softConstraintVals.size(), basisCyclesMatrix.cols());
    withAllConstraintMatrix.setFromTriplets(softConstraintMatrixTriplets.begin(), softConstraintMatrixTriplets.end());

    Eigen::VectorXd actualCurvatureAllConstraints = withAllConstraintMatrix * angleJumpsPerEdge;
    Eigen::VectorXd diff = withAllConstraintMatrix * edgeJumps;

    cycleCurvature.conservativeResize(basisCyclesMatrix.rows() + hardConstraintVals.size() + softConstraintVals.size());
    cycleIndices.conservativeResize(basisCyclesMatrix.rows() + hardConstraintVals.size() + softConstraintVals.size());

    for (int i = 0; i < softConstraintVals.size(); i++) {
      int id = basisCyclesMatrix.rows() + hardConstraintVals.size() + i;
      double constraint = softConstraintVals[i];
      double naturalDiff = diff(id) + actualCurvatureAllConstraints(id);

      if (constraint - naturalDiff > igl::PI) {
        while (constraint - naturalDiff > igl::PI) {
          constraint -= 2 * igl::PI;
        }
      } else if (naturalDiff - constraint > igl::PI) {
        while (naturalDiff - constraint > igl::PI) {
          constraint += 2 * igl::PI;
        }
      }

      // if (i == 6) {
      //   constraint += 2 * igl::PI;
      // }
      // if (i == 5) {
      //   constraint += 2 * igl::PI;
      // }
      // if (i == 4) {
      //   constraint -= 2 * igl::PI;
      // }

      std::cout << constraint << ", " << naturalDiff << ", " << softConstraintVals[i] << std::endl;

      cycleCurvature(id) = actualCurvatureAllConstraints(id);
      cycleIndices(id) = constraint / (2.0 * igl::PI / (double)N);
    }

    matrix = withAllConstraintMatrix;
    rp::index_prescription(V, F, EV, EF, innerEdges, withAllConstraintMatrix, cycleCurvature, cycleIndices, solver, N, weightMatrix, edgeJumps, linfError);
  }
}

