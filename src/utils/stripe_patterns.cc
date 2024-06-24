#include "./stripe_patterns.h"

#include "./vector_utils.h"

#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <iostream>

namespace rp {
namespace stripe_patterns {
namespace {
  double EPSILON = 1e-7;
}

double norm(std::complex<double> &cplx) {
  return std::sqrt(cplx.real() * cplx.real() + cplx.imag() * cplx.imag());
}

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::VectorXd &vectorAnglesPerVertex,
             const Eigen::VectorXd &rescalings, const double &velocity,
             Eigen::VectorXd &parameterization) {
  Eigen::VectorXd complexParameterization(2 * V.rows());
  
  Eigen::SparseMatrix<double> energyMatrix(2 * V.rows(), 2 * V.rows());
  std::vector<Eigen::Triplet<double>> energyMatrixTriplets = {};

  Eigen::VectorXd cotanWeights(EV.rows());
  cotanWeights.setZero();
  Eigen::VectorXcd complexParameterizationVals(V.rows());
  Eigen::VectorXcd complexEdgeJumps(EV.rows());

  // cotan weightを計算
  for (int i = 0; i < EV.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == 0) {
      continue;
    }

    Eigen::Vector2d cotans(2);
    cotans.setZero();
    for (int j = 0; j < 2; j++) {
      int fid = EF(i, j);
      if (fid == -1) {
        continue;
      }
      for (int k = 0; k < 3; k++) {
        int vid = F(fid, k);
        if (vid != EV(i, 0) && vid != EV(i, 1)) {
          double angle = faceAngles(fid, k);
          if (std::abs(angle) > EPSILON * 10000) {
            double cot = std::cos(angle) / std::sin(angle);
            cotans(j) = std::cos(angle) / std::sin(angle);
          }
          break;
        }
      }
    }

    cotanWeights(i) = std::min(std::max((cotans(0) + cotans(1)) / 2, 0.001), 10.);
  }

  // 対角要素以外
  for (int i = 0; i < EV.rows(); i++) {
    // 必ずv_i < v_jとなるようにする
    int v_i = EV(i, 0) < EV(i, 1) ? EV(i, 0) : EV(i, 1);
    int v_j = EV(i, 0) < EV(i, 1) ? EV(i, 1) : EV(i, 0);
    double length = (V.row(v_i) - V.row(v_j)).norm();

    Eigen::Vector2d theta;
    theta.setZero();
    for (int j = 0; j < 2; j++) {
      int vid = EV(i, j);
      bool visited = false;
      for (int k = 0; k < VE[vid].size(); k++) {
        if (VE[vid][k] == i) {
          theta(j) = VEAngles[vid][k];
          visited = true;
          // v_i->v_jの向きが正なので、もしvidがv_jならπ回転させる
          if (vid == v_j) {
            theta(j) += igl::PI;
          }
          break;
        }
      }
      assert(visited);
    }

    double arg0 = velocity * rescalings(EV(i, 0)) * std::cos(vectorAnglesPerVertex(EV(i, 0)) - theta(0));
    double arg1 = velocity * rescalings(EV(i, 1)) * std::cos(vectorAnglesPerVertex(EV(i, 1)) - theta(1));

    // edgeJumpが90度を超える場合は1を超える場合は共役を取らなければならないらしいが、そんなことは多分ない
    double angle = length * (arg0 + arg1) / 2;

    complexEdgeJumps(i) = std::exp(std::complex<double>(0, 1) * angle);

    Eigen::Matrix2d mat;
    mat << std::cos(angle), std::sin(angle), -std::sin(angle), std::cos(angle);
    mat = -cotanWeights(i) * mat;
    Eigen::Matrix2d tranMat = mat.transpose();

    for (int k = 0; k < 2; k++) {
      for (int l = 0; l < 2; l++) {
        energyMatrixTriplets.emplace_back(Eigen::Triplet<double>(2 * v_i + k, 2 * v_j + l, mat(k, l)));
        energyMatrixTriplets.emplace_back(Eigen::Triplet<double>(2 * v_j + k, 2 * v_i + l, tranMat(k, l)));
      }
    }
  }

  for (int i = 0; i < EV.rows(); i++) {
    for (int k = 0; k < 2; k++) {
      int vid = EV(i, k);
      for (int l = 0; l < 2; l++)
        energyMatrixTriplets.emplace_back(Eigen::Triplet<double>(2 * vid + l, 2 * vid + l, cotanWeights(i)));
    }
  }

  energyMatrix.setFromTriplets(energyMatrixTriplets.begin(), energyMatrixTriplets.end());

  Eigen::SparseMatrix<double> weightMatrix(2 * V.rows(), 2 * V.rows());
  std::vector<Eigen::Triplet<double>> weightMatrixTriplets = {};
  double areaWeight = std::sqrt(doubleAreas.sum() / 6);
  // barycentric weightなので各三角形の1/3の面積を押し付ける
  for (int i = 0; i < F.rows(); i++) {
    double area = doubleAreas(i) / 6;
    for (int j = 0; j < 3; j++) {
      int vid = F(i, j);
      for (int k = 0; k < 2; k++) {
        weightMatrixTriplets.emplace_back(Eigen::Triplet<double>(2 * vid + k, 2 * vid + k, area / areaWeight));
      }
    }
  }

  weightMatrix.setFromTriplets(weightMatrixTriplets.begin(), weightMatrixTriplets.end());

  std::cout << cotanWeights.maxCoeff() << ", " << cotanWeights.minCoeff() << std::endl;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
  ldltSolver.compute(energyMatrix);

  assert(ldltSolver.info() == Eigen::Success);

  complexParameterization = Eigen::VectorXd::Random(2 * V.rows(), 1);

  for (int i = 0; i < 30; i++) {
    complexParameterization = ldltSolver.solve(weightMatrix * complexParameterization);
    double norm = complexParameterization.transpose() * weightMatrix * complexParameterization;
    complexParameterization = complexParameterization / std::sqrt(norm);
  }

  std::cout << "energy: " << (complexParameterization.transpose() * energyMatrix * complexParameterization) << ", xTBx: " << complexParameterization.transpose() * weightMatrix * complexParameterization << ", B sum: " << weightMatrix.sum() << std::endl;

  parameterization.resize(V.rows());

  Eigen::VectorXd errors(EV.rows());
  errors.setZero();

  for (int i = 0; i < V.rows(); i++) {
    double angle = std::atan2(complexParameterization(2 * i + 1), complexParameterization(2 * i));
    complexParameterizationVals(i) = std::complex<double>((complexParameterization(2 * i), complexParameterization(2 * i + 1)));
    parameterization(i) = round_pi(angle);
  }

  std::cout << "actual energy: " << complexParameterization.transpose() * energyMatrix * complexParameterization << std::endl;
};
}
}
