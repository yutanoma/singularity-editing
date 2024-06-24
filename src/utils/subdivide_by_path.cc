#include "./subdivide_by_path.h"

#include <igl/is_edge_manifold.h>
#include <igl/boundary_loop.h>

namespace rp {
namespace subdivide_by_path {
namespace {
  double EPSILON = 1e-7;
}

void process(const Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, std::vector<std::vector<std::array<int, 2>>> &pathOnFaces) {
  std::vector<std::array<int, 3>> newFaces = {};
  std::cout << "subdivision" << std::endl;

  for (int i = 0; i < pathOnFaces.size(); i++) {
    auto paths = pathOnFaces[i];
    auto pathsNum = paths.size();
    int fid = i;

    // fixme: 現状はpathは一つの三角形メッシュあたり1つしか認めてない
    assert(paths.size() <= 1);

    for (int j = 0; j < paths.size(); j++) {
      int e0id = paths[j][0];
      int e1id = paths[j][1];
      Eigen::Vector3d edge0 = V.row(e0id).transpose();
      Eigen::Vector3d edge1 = V.row(e1id).transpose();

      bool divided = false;

      // どのedgeに乗っているかを判定し、各edgeの共有点とその点を結ぶ
      for (int k = 0; k < 3; k++) {
        int vid = F(fid, k);
        int v0id = F(fid, (k + 1) % 3);
        int v1id = F(fid, (k + 2) % 3);

        // e0またはe1がv, v0, v1のいずれかと一致していればその時点で終了
        if (k == 0) {
          int correspondingVid = -1;
          for (int l = 0; l < 2; l++) {
            int edgeVid = paths[j][l];
            int otherEdgeVid = paths[j][(l + 1) % 2];
            if (vid == edgeVid) {
              F.row(fid) << vid, v0id, otherEdgeVid;
              std::array<int, 3> nf = {vid, otherEdgeVid, v1id};
              newFaces.emplace_back(nf);
              correspondingVid = vid;
              break;
            } else if (v0id == edgeVid) {
              F.row(fid) << v0id, v1id, otherEdgeVid;
              std::array<int, 3> nf = {v0id, otherEdgeVid, vid};
              newFaces.emplace_back(nf);
              correspondingVid = v0id;
              break;
            } else if (v1id == edgeVid) {
              F.row(fid) << v1id, vid, otherEdgeVid;
              std::array<int, 3> nf = {v1id, otherEdgeVid, v0id};
              newFaces.emplace_back(nf);
              correspondingVid = v1id;
              break;
            }
          }

          if (correspondingVid != -1) {
            divided = true;
            break;
          }
        }

        Eigen::Vector3d v = V.row(vid).transpose();
        Eigen::Vector3d otherV0 = V.row(v0id).transpose();
        Eigen::Vector3d otherV1 = V.row(v1id).transpose();

        Eigen::Vector3d v2edge0 = edge0 - v;
        Eigen::Vector3d v2edge1 = edge1 - v;
        Eigen::Vector3d v02edge0 = edge0 - otherV0;
        Eigen::Vector3d v12edge0 = edge0 - otherV1;
        Eigen::Vector3d v02edge1 = edge1 - otherV0;
        Eigen::Vector3d v12edge1 = edge1 - otherV1;
        Eigen::Vector3d v2v0 = otherV0 - v;
        Eigen::Vector3d v2v1 = otherV1 - v;

        if (
          std::abs(v2edge0.norm() + v02edge0.norm() - v2v0.norm()) < EPSILON
          && std::abs(v2edge1.norm() + v12edge1.norm() - v2v1.norm()) < EPSILON
        ) {
          // otherVid0とvidの間にe0があり、かつotherVid1とvidの間にe1があればOK
          // (v, e0, e1)(e0, v0, v1)(e1, e0, v1)の三角形が新たにできる
          F.row(fid) << vid, e0id, e1id;
          std::array<int, 3> nf1 = {e0id, v0id, v1id},
                             nf2 = {e1id, e0id, v1id};
          newFaces.emplace_back(nf1);
          newFaces.emplace_back(nf2);
          divided = true;
          // std::cout << fid << ": type1" << std::endl;
          break;
        }
        else if (
          std::abs(v2edge1.norm() + v02edge1.norm() - v2v0.norm()) < EPSILON
          && std::abs(v2edge0.norm() + v12edge0.norm() - v2v1.norm()) < EPSILON
        ) {
          // otherVid0とvidの間にe1があり、かつotherVid1とvidの間にe0があればOK
          // (e0, v, e1)(v0, e0, e1)(v0, v1, e0)の三角形が新たにできる
          F.row(fid) << e0id, vid, e1id;
          std::array<int, 3> nf1 = {v0id, e0id, e1id},
                             nf2 = {v0id, v1id, e0id};
          newFaces.emplace_back(nf1);
          newFaces.emplace_back(nf2);
          divided = true;
          // std::cout << fid << ": type2";
          break;
        }
      }

      assert(divided);
    }
  }

  int originalFacesNum = F.rows();
  F.conservativeResize(originalFacesNum + newFaces.size(), 3);
  for (int i = 0; i < newFaces.size(); i++) {
    for (int j = 0; j < 3; j++) {
      F(originalFacesNum + i, j) = newFaces[i][j];
    }
  }

  std::cout << "subdivision finished" << std::endl;

  // std::vector<std::vector<int>> loop1;
  // igl::boundary_loop(F, loop1);
  // std::cout << loop1.size() << std::endl;

  // std::cout << igl::is_edge_manifold(F) << std::endl;
};
}
}
