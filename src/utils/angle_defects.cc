#include "angle_defects.h"

#include "./vector_utils.h"

#include <igl/PI.h>
#include <iostream>

namespace rp {
void angle_defects(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                   const Eigen::VectorXd &anglesSumPerVertex,
                   const Eigen::MatrixX3d &faceAngles,
                   const Eigen::SparseMatrix<double> &basisLoops,
                   const std::vector<std::vector<int>> &boundaryLoops,
                   const Eigen::VectorXd &angleJumps,
                   Eigen::VectorXd &angleDefects) {
  angleDefects.resize(basisLoops.rows());

  angleDefects = basisLoops * angleJumps;

  for (int i = 0; i < F.rows(); i++) {
    double sum = 0;
    for (int j = 0; j < 3; j++) {
      double angleSum = anglesSumPerVertex(F(i, j));
      double angle = faceAngles(i, j);
      sum += angle * 2 * igl::PI / angleSum;
    }
    double oldAngle = angleDefects(i);
    angleDefects(i) = sum - igl::PI;
    // if (i < 10) {
    //   for (int j = 0; j < basisLoops.cols(); j++) {
    //     if (basisLoops.coeff(i, j) != 0)
    //       std::cout << "(" << basisLoops.coeff(i, j) << ", " << i << ", " << j << ") ";
    //   }
    //   std::cout << std::endl << oldAngle << " " << angleDefects(i) << std::endl;
    // }
  }

  Eigen::VectorXd vertexAnglesSum(V.rows());
  vertexAnglesSum.setZero();

  for (int i = 0; i < faceAngles.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      vertexAnglesSum(F(i, j)) += faceAngles(i, j);
    }
  }

  for (int i = 0; i < boundaryLoops.size(); i++) {
    double sum = 0;
    for (int j = 0; j < boundaryLoops[i].size(); j++) {
      int vid = boundaryLoops[i][j];
      double anglesSum = vertexAnglesSum(vid);
      sum += anglesSum - igl::PI;
    }

    int idx = F.rows() + i;
    double oldVal = angleDefects(idx);
    angleDefects(idx) = 2 * igl::PI - sum;
    // std::cout << oldVal << ", " << angleDefects(idx) << std::endl;
  }

  for (int i = 0; i < (basisLoops.rows() - F.rows() - boundaryLoops.size()); i++) {
    int idx = F.rows() + boundaryLoops.size() + i;
    angleDefects(idx) = rp::round_pi(angleDefects(idx));
  }
}
}
