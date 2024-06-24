#include "./curl_correction.h"

#include "./vector_utils.h"

#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <iostream>

namespace rp {
namespace curl_correction {

// faceAngles: faceのある頂点の角度がいくつかどうか
// VEAngles: あるvertexに隣接するedgeへの2πに丸めた角度がいくつかどうか
// vectorAnglesPerVertex: ある角度でベクトル場の回転角がいくつか。inner vertexでは合計が2πに丸められている
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixX3d &B1,
             const Eigen::MatrixX3d &B2, const Eigen::MatrixX3d &B3,
             const std::vector<std::vector<int>> &VE, const std::vector<std::vector<double>> &VEAngles,
             const Eigen::VectorXd &doubleAreas,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             const Eigen::MatrixX3d &faceAngles,
             const Eigen::VectorXd &vectorAnglesPerVertex,
             const Eigen::VectorXd &anglesSumPerVertex,
            Eigen::VectorXd &rescalings) {
  Eigen::VectorXd resultVec(2 * F.rows()), weightVec(2 * F.rows());
  // curl correctionの最小二乗法の行列。F.rows()*2行ある。
  Eigen::SparseMatrix<double> weightMatrix(2 * F.rows(), V.rows() - 1);
  std::vector<Eigen::Triplet<double>> weightMatrixTriplets = {};

  std::cout << "curl correction start" << std::endl;

  Eigen::MatrixX3d gammaVals(F.rows(), 3);

  for (int i = 0; i < F.rows(); i++) {
    Eigen::Matrix<double, 2, 3> jacobian(2, 3);
    Eigen::Vector3d gammas;

    Eigen::Matrix<double, 3, 2> xyList(3, 2);

    for (int j = 0; j < 3; j++) {
      // j = 0の時は座標は(.0, .0)
      if (j == 0) {
        xyList.row(j) << .0, .0;
      } else {
        Eigen::Vector3d vect = (V.row(F(i, j)) - V.row(F(i, 0))).transpose();
        double x = vect.dot(B1.row(F(i, j)).transpose());
        double y = vect.dot(B2.row(F(i, j)).transpose());
        xyList.row(j) << x, y;
      }
    }

    // jacobianとgammasに実際の値を入れる
    jacobian << xyList(1, 1) - xyList(2, 1), xyList(2, 1) - xyList(0, 1), xyList(0, 1) - xyList(1, 1),
                xyList(2, 0) - xyList(1, 0), xyList(0, 0) - xyList(2, 0), xyList(1, 0) - xyList(0, 0);

    Eigen::Vector3d triangleAngles = faceAngles.row(i).transpose();

    for (int j = 0; j < 3; j++) {
      // boundary上ではs_iは1。[Knoppel+, Globally Optimal Direction Fields, 2013]参照
      double angle = vectorAnglesPerVertex(F(i, j));

      double theta = 0;
      // F(i, j)->F(i, j+1)の辺を探す
      for (int k = 0; k < VE[F(i, j)].size(); k++) {
        const int edgeId = VE[F(i, j)][k];
        const int vid = F(i, j);
        const int nextVid = F(i, (j + 1) % 3);
        if (
          (EV(edgeId, 0) == vid && EV(edgeId, 1) == nextVid) ||
          (EV(edgeId, 1) == vid && EV(edgeId, 0) == nextVid)
        ) {
          double angle2Edge = VEAngles[F(i, j)][k];
          if (j == 0) {
            theta = angle - angle2Edge;
          } else if (j == 1) {
            theta = angle - (angle2Edge + triangleAngles(1) - igl::PI);
          } else if (j == 2) {
            theta = angle - (angle2Edge + triangleAngles(1) + triangleAngles(2));
          }
          break;
        }
      }

      double gamma = theta;
      gammas(j) = round_pi(gamma);
    }

    // 最大と最小がπ以上離れていたら明らかに何かがおかしいので、サンプリング範囲を変更する
    if (gammas.maxCoeff() - gammas.minCoeff() > igl::PI) {
      for (int i = 0; i < 3; i++) {
        if (gammas(i) < 0) {
          gammas(i) = gammas(i) + 2 * igl::PI;
        }
      }
    }

    // if (F(i, 0) == 8083 || F(i, 1) == 8083 || F(i, 2) == 8083) {
    //   std::cout << gammas.transpose() << ", " << doubleAreas(i) << std::endl;
    // }

    gammaVals.row(i) << gammas.transpose();

    jacobian = jacobian * doubleAreas(i) / 4;

    Eigen::Matrix2d rot;
    rot << 0, -1, 1, 0;
    Eigen::Vector2d b = rot * jacobian * gammas;

    double weight = std::sqrt(doubleAreas(i) / 2);

    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 3; k++) {
        if (F(i, k) == 0) {
          // vidが0の場合はlog(ν)は0になる
          continue;
        }
        Eigen::Triplet<double> triplet(2 * i + j , F(i, k) - 1, jacobian(j, k) * weight);
        weightMatrixTriplets.emplace_back(triplet);
      }
      resultVec(2 * i + j) = b(j) * weight;
      weightVec(2 * i + j) = weight;
    }
  }

  weightMatrix.setFromTriplets(weightMatrixTriplets.begin(), weightMatrixTriplets.end());

  std::cout << "solve curl correction" << std::endl;

  // solve the least square solution
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
  ldltSolver.compute(weightMatrix.transpose() * weightMatrix);

  assert(ldltSolver.info() == Eigen::Success);

  Eigen::VectorXd rescalingLogs(V.rows());
  Eigen::VectorXd _rescalingLogs = ldltSolver.solve(weightMatrix.transpose() * resultVec);

  // std::cout << _rescalingLogs << std::endl;

  std::cout << "solved, linfError: " << (weightMatrix.transpose() * weightMatrix * _rescalingLogs - weightMatrix.transpose() * resultVec).norm() << ", energy: " << (weightMatrix * _rescalingLogs - resultVec).norm() << std::endl;

  rescalingLogs << .0, _rescalingLogs;

  rescalings.resize(rescalingLogs.rows());

  for (int i = 0; i < rescalingLogs.rows(); i++) {
    rescalings(i) = std::exp(rescalingLogs(i));
  }

  // for (int i = 0; i < F.rows(); i++) {
  //   std::cout << "face[" << i << "]: " << gammaVals.row(i) << ", " << doubleAreas(i) << ", [";
  //   for (int j = 0; j < 3; j++) {
  //     int vid = F(i, j);
  //     std::cout << "[" << rescalings(vid) << ", " << rescalingLogs(vid) << "], ";
  //   }
  //   std::cout << "]" << std::endl;
  // }

  double mean = rescalings.mean();

  std::cout << rescalings.maxCoeff() << ", " << rescalings.minCoeff() << ", " << mean << std::endl;
};
}
}

