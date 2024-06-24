#include "hodge_star.h"

#include <Eigen/Dense>

namespace rp {
void hodge_star(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV, const Eigen::VectorXi &innerEdges, Eigen::SparseMatrix<double> &hodgeStar) {
  hodgeStar.resize(innerEdges.rows(), innerEdges.rows());
  hodgeStar.setZero();
  std::vector<Eigen::Triplet<double>> hodgeStarTripletVec;
  for (int i = 0; i < innerEdges.rows(); i++) {
    double cotangentWeight = 0;
    // cotangent weightを計算する
    for (int j = 0; j < 2; j++) {
      int fid = EF(innerEdges(i), j);
      if (fid == -1)
        continue;
      // このfidの内、edgeを構成していない点を探す
      int otherVid = -1, otherVidIdx = -1;
      for (int k = 0; k < 3; k++) {
        if (F(fid, k) != EV(innerEdges(i), 0) && F(fid, k) != EV(innerEdges(i), 1)) {
          otherVid = F(fid, k);
          otherVidIdx = k;
        }
      }
      assert(otherVid != -1);

      // otherVidの角度のcotangentを計算する
      Eigen::Vector3d vec1 = V.row(F(fid, (otherVidIdx + 1) % 3)).transpose() - V.row(otherVid).transpose();
      Eigen::Vector3d vec2 = V.row(F(fid, (otherVidIdx + 2) % 3)).transpose() - V.row(otherVid).transpose();
      double angle = rp::get_angle(vec1, vec2, vec1.cross(vec2));
      cotangentWeight += std::cos(angle) / std::sin(angle);
    }
    hodgeStarTripletVec.emplace_back(Eigen::Triplet<double>(i, i, std::sqrt(std::max(cotangentWeight / 2, 0.001))));
    // hodgeStarTripletVec.emplace_back(Eigen::Triplet<double>(i, i, 1.0));
  }
  hodgeStar.setFromTriplets(hodgeStarTripletVec.begin(), hodgeStarTripletVec.end());
}
}