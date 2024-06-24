#include "./vector_utils.h"

#include <Eigen/LU>
#include <Eigen/Dense>

namespace rp {
double get_angle(const Eigen::Vector3d& basis_vec,
                 const Eigen::Vector3d& angled_vec,
                 const Eigen::Vector3d& normal) {
  double angle =
      atan2(normal.normalized().dot(basis_vec.normalized().cross(angled_vec.normalized())), angled_vec.normalized().dot(basis_vec.normalized()));

  assert(angle >= -igl::PI && angle <= igl::PI);

  return angle;
}

Eigen::Vector3d face_normal(const Eigen::MatrixX3d& V,
                            const Eigen::MatrixX3i& F, const int face_id) {
  Eigen::Vector3d v0 =
      V.row(F(face_id, 1)).transpose() - V.row(F(face_id, 0)).transpose();
  Eigen::Vector3d v1 =
      V.row(F(face_id, 2)).transpose() - V.row(F(face_id, 0)).transpose();

  return v0.cross(v1).normalized();
}

// 2次元ベクトルxを、v_aとv_bの線形和で表す
// result(0)がv_a成分、result(1)がv_b成分
void vector_decomposition(Eigen::Vector2d &x, Eigen::Vector2d &v_a,
                          Eigen::Vector2d &v_b,
                          Eigen::Vector2d &result) {
  Eigen::Matrix2d matrix;
  matrix << v_a(0), v_b(0), v_a(1), v_b(1);
  result = matrix.inverse() * x;
};

// 3次元ベクトルxを、v_aとv_bの線形和で表す
// v_aとv_bの張る平面上にない場合にはなぞの値が出てくる。
void vector_decomposition_3d(Eigen::Vector3d &x, Eigen::Vector3d &v_a,
                             Eigen::Vector3d &v_b,
                             Eigen::Vector2d &result) {
  Eigen::Matrix2d matrix;
  matrix << v_a(0), v_b(0), v_a(1), v_b(1);
  Eigen::Vector2d _x;
  _x << x(0), x(1);
  result = matrix.inverse() * _x;
};

void get_edge_from_faces(const int& prevFid, const int& postFid,
                         const Eigen::MatrixXi& FE, const Eigen::MatrixXi& EF,
                         int& edge, int& edgeSign) {
  int edgeId = -1;

  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      if (FE(prevFid, j) == FE(postFid, k)) {
        edgeId = FE(prevFid, j);
      }
    }
  }

  assert(edgeId != -1 &&
         (EF(edgeId, 0) == prevFid || EF(edgeId, 0) == postFid) &&
         (EF(edgeId, 1) == prevFid || EF(edgeId, 1) == postFid));

  edge = edgeId;

  // 0->1の方向が正
  if (EF(edgeId, 0) == prevFid) {
    edgeSign = 1;
  } else if (EF(edgeId, 1) == prevFid) {
    edgeSign = -1;
  } else {
    assert("invalid edge from faces");
  }
};

double round_pi(const double &angle) {
  double ad = angle;
  while (ad <= - igl::PI || ad > igl::PI) {
    if (ad > igl::PI) {
      ad = ad - 2 * igl::PI;
    } else {
      ad = ad + 2 * igl::PI;
    }
  }

  return ad;
}

void create_tree_from_tEf(
    const Eigen::VectorXi& dualTreeFathers, const Eigen::MatrixXi& EV,
    std::vector<std::vector<int>>& treeNextNodeList) {
  treeNextNodeList.resize(dualTreeFathers.rows(), {});

  for (int i = 0; i < dualTreeFathers.rows(); i++) {
    if (dualTreeFathers(i) != -1) {
      int other_vid = -1;

      for (int j = 0; j < 2; j++) {
        if (EV(dualTreeFathers(i), j) != i) {
          other_vid = EV(dualTreeFathers(i), j);
          break;
        }
      }

      assert(other_vid != -1);

      treeNextNodeList[other_vid].emplace_back(i);
    }
  }
}
}  // namespace rp
