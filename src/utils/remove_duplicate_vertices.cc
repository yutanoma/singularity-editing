#include "remove_duplicate_vertices.h"

namespace rp {
void remove_duplicate_vertices(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F) {
  Eigen::VectorXi J;
  remove_duplicate_vertices(V, F, J, 0.0000001);
}

void remove_duplicate_vertices(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, Eigen::VectorXi &SVJ, const double &val) {
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SVI;
  // 重複頂点の削除
  igl::remove_duplicate_vertices(V, val, SV, SVI, SVJ);
  // SVに削除後の頂点
  V = SV;
  // SVJに元の頂点番号→削除後の頂点番号が格納されているので、
  // インデックスリストを更新
  // ループを回さず一括処理する方法があるはず？
  Eigen::VectorXi isDegenerate(F.rows());
  isDegenerate.setZero();
  Eigen::MatrixX3i NF(F);
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < F.cols(); j++) {
      int temp = NF(i, j);
      NF(i, j) = SVJ(temp);
    }
  }

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      if (NF(i, j) == NF(i, (j + 1) % 3)) {
        isDegenerate(i) = 1;
      }
    }
  }

  F.resize(NF.rows() - isDegenerate.sum(), 3);
  int count = 0;
  for (int i = 0; i < NF.rows(); i++) {
    if (!isDegenerate(i)) {
      F.row(count) = NF.row(i);
      count++;
    }
  }
};
}  // namespace rp