#include "cut_mesh_with_path.h"

#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>

#include "./vector_utils.h"

namespace rp {
namespace cut_mesh_with_path {
namespace {
double EPSILON = 1e-7;

enum PointType
{
  point,
  edge
};
}

void computeVF(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V, std::vector<std::vector<int>> &VF) {
  VF.clear();
  VF.resize(V.rows(), {});

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < F.cols(); j++) {
      VF[F(i, j)].emplace_back(i);
    }
  }
}

// pathsは必ず元のメッシュの上に乗っているものとする。乗っていない場合には正常に動作しない
// またループには対応していない。開曲線のみ
void process(const Eigen::MatrixX3i &F, Eigen::MatrixX3d &V, const std::vector<Eigen::MatrixX3d> &paths, std::vector<std::vector<std::array<int, 2>>> &vidsPerFace, std::vector<std::vector<int>> &pathVids) {
  if (paths.size() == 0) {
    return;
  }

  Eigen::MatrixXi EV, FE, EF;
  igl::edge_topology(V, F, EV, FE, EF);

  std::vector<std::vector<int>> VF = {};

  computeVF(F, V, VF);

  // 今追加されているvidsの中で最も大きなもの
  int frontierVid = V.rows() - 1;
  // 新しいvertices
  std::vector<Eigen::Vector3d> newVertices = {};

  vidsPerFace.resize(F.rows(), {});

  pathVids.resize(paths.size(), {});

  for (int i = 0; i < paths.size(); i++) {
    Eigen::Vector3d firstPoint = paths[i].row(0).transpose();
    auto path = paths[i];

    // 最も近い点を探す
    int nearestPointVid = -1;
    double distance = 99999;
    for (int j = 0; j < V.rows(); j++) {
      Eigen::Vector3d p = V.row(j).transpose();
      auto d = (p - firstPoint).norm();
      if (distance > d) {
        distance = d;
        nearestPointVid = j;
      }

      if (distance < EPSILON) {
        break;
      }
    }

    assert(nearestPointVid != -1);

    int prevVid = nearestPointVid;
    int prevEdgeId = -1;
    PointType prevPointType = PointType::point;

    int currentFid = -1;

    pathVids[i].emplace_back(nearestPointVid);

    std::vector<int> pathFids = {};

    // faceIdからパスの方向に伸ばしていく
    for (int j = 1; j < path.rows(); j++) {
      Eigen::Vector3d prevPoint = path.row(j - 1).transpose();
      Eigen::Vector3d nextPoint = path.row(j).transpose();
      double prevNextDistance = (prevPoint - nextPoint).norm();
      double coeffSum = 0;

      while (true) {
        // 次の点を探す。
        if (prevPointType == PointType::point) {
          // まず次の平面を探す
          auto adjacentFids = VF[prevVid];

          for (int i = 0; i < adjacentFids.size(); i++) {
            int fid = adjacentFids[i];

            Eigen::Vector3d va = V.row(F(fid, 0)).transpose();
            Eigen::Vector3d vb = V.row(F(fid, 1)).transpose();
            Eigen::Vector3d vc = V.row(F(fid, 2)).transpose();

            Eigen::Vector3d ab = vb - va;
            Eigen::Vector3d ac = vc - va;

            Eigen::Vector3d ad = prevPoint - va;
            Eigen::Vector3d ae = nextPoint - va;

            // ab, ac, ad, aeがすべて同一平面上にある
            Eigen::Matrix3d abcd;
            abcd << ab(0), ab(1), ab(2),
                    ac(0), ac(1), ac(2),
                    ad(0), ad(1), ad(2);
            
            Eigen::Matrix3d abce;
            abcd << ab(0), ab(1), ab(2),
                    ac(0), ac(1), ac(2),
                    ae(0), ae(1), ae(2);

            // 同一平面上にない場合はスルー
            if (std::abs(abcd.determinant()) > EPSILON || std::abs(abce.determinant()) > EPSILON) {
              continue;
            }

            // 係数が両方0以上でなければスルー
            // 前の点がpointだった場合
            // vector_decompositionをして、残りの2点の作る辺の上の点を求める
            int vaId, vbId;
            for (int i = 0; i < 3; i++) {
              if (F(fid, i) == prevVid) {
                vaId = F(fid, (i + 1) % 3);
                vbId = F(fid, (i + 2) % 3);
                va = V.row(vaId).transpose() - prevPoint;
                vb = V.row(vbId).transpose() - prevPoint;
              }
            }

            Eigen::Vector2d result;
            Eigen::Vector3d vec = nextPoint - prevPoint;
            rp::vector_decomposition_3d(vec, va, vb, result);

            if (result(0) < -EPSILON || result(1) < -EPSILON) {
              continue;
            }

            currentFid = fid;
            break;
          }

          assert(currentFid != -1);

          // 前の点がpointだった場合
          // vector_decompositionをして、残りの2点の作る辺の上の点を求める
          Eigen::Vector3d va, vb;
          int vaId, vbId;
          for (int i = 0; i < 3; i++) {
            if (F(currentFid, i) == prevVid) {
              vaId = F(currentFid, (i + 1) % 3);
              vbId = F(currentFid, (i + 2) % 3);
              va = V.row(vaId).transpose() - prevPoint;
              vb = V.row(vbId).transpose() - prevPoint;
            }
          }

          Eigen::Vector2d result;
          Eigen::Vector3d vec = nextPoint - prevPoint;
          rp::vector_decomposition_3d(vec, va, vb, result);

          coeffSum += (result(0) * va + result(1) * vb).norm() / prevNextDistance;

          if (std::abs(result(0)) < EPSILON ^ std::abs(result(1)) < EPSILON) {
            // 次も点
            if (std::abs(result(0)) < EPSILON) {
              // result(1)=1、つまりvbが次の点
              prevVid = vbId;
            } else {
              // result(0)=1、つまりvaが次の点
              prevVid = vaId;
            }
            prevPointType = PointType::point;
          } else {
            // 次は辺
            // 新しい点を追加する処理
            Eigen::Vector3d newPoint = prevPoint + result(0) * va + result(1) * vb;
            newVertices.emplace_back(newPoint);
            frontierVid += 1;
            std::array<int, 2> facePath = {prevVid, frontierVid};
            vidsPerFace[currentFid].emplace_back(facePath);

            // newPointの乗るedgeのid
            for (int i = 0; i < 3; i++) {
              int edgeId = FE(currentFid, i);
              if ((vaId == EV(edgeId, 0) && vbId == EV(edgeId, 1)) || (vaId == EV(edgeId, 1) && vbId == EV(edgeId, 0))) {
                prevEdgeId = edgeId;
                break;
              }
            }

            prevVid = frontierVid;
            prevPointType = PointType::edge;
          }

        } else if (prevPointType == PointType::edge) {
          // 前の点がedgeだった場合
          assert(prevEdgeId != -1);

          currentFid = EF(prevEdgeId, 0) == currentFid ? EF(prevEdgeId, 1) : EF(prevEdgeId, 0);
          Eigen::Vector2i edgeVids = EV.row(prevEdgeId).transpose();
          int otherVid = -1;

          for (int i = 0; i < 3; i++) {
            if (F(currentFid, i) != edgeVids(0) && F(currentFid, i) != edgeVids(1)) {
              otherVid = F(currentFid, i);
              break;
            }
          }
          assert(otherVid != -1);

          bool foundFlag = false;
          for (int i = 0; i < 2; i++) {
            int edgeVid = edgeVids(i);
            // vaをedgeVidとの間、vbをotherVidとの間のベクトルとし、vector_decompositionをする
            // edgeの最新行が常に直掩のedgeの座標になっているはず
            Eigen::Vector3d prevV = newVertices[newVertices.size() - 1];
            Eigen::Vector3d va = V.row(edgeVid).transpose() - prevV;
            Eigen::Vector3d vb = V.row(otherVid).transpose() - prevV;

            Eigen::Vector3d vec = nextPoint - prevV;

            Eigen::Vector2d result;
            rp::vector_decomposition_3d(vec, va, vb, result);

            coeffSum += (result(0) * va + result(1) * vb).norm() / prevNextDistance;

            if (std::abs(result(0)) < EPSILON ^ std::abs(result(1)) < EPSILON) {
              // 次は点
              int evid = prevVid;

              if (std::abs(result(0)) < EPSILON) {
                // result(1)=1、つまりvbが次の点
                prevVid = otherVid;
              } else {
                // result(0)=1、つまりvaが次の点
                prevVid = edgeVid;
              }

              std::array<int, 2> facePath = {evid, prevVid};
              vidsPerFace[currentFid].emplace_back(facePath);

              prevPointType = PointType::point;

              foundFlag = true;
              break;
            } else if (result(0) > -EPSILON && result(1) > EPSILON) {
              // 次も辺
              // 新しい点を追加する処理
              Eigen::Vector3d newPoint = prevPoint + result(0) * va + result(1) * vb;
              newVertices.emplace_back(newPoint);
              frontierVid += 1;
              std::array<int, 2> facePath = {prevVid, frontierVid};
              vidsPerFace[currentFid].emplace_back(facePath);

              // newPointの乗るedgeのid
              for (int i = 0; i < 3; i++) {
                int edgeId = FE(currentFid, i);
                if ((otherVid == EV(edgeId, 0) && edgeVid == EV(edgeId, 1)) || (edgeVid == EV(edgeId, 0) && otherVid == EV(edgeId, 1))) {
                  prevEdgeId = edgeId;
                  break;
                }
              }

              prevVid = frontierVid;
              prevPointType = PointType::edge;

              foundFlag = true;
              break;
            }
          }

          assert(foundFlag);
        }

        pathVids[i].emplace_back(prevVid);

        pathFids.emplace_back(currentFid);

        if (coeffSum - 1 > - EPSILON) {
          // 係数を足して1なら、それは目的の点に到達したということなので次へ行く。
          break;
        }

      }
    }

    for (int j = 0; j < pathVids[i].size(); j++) {
      std::cout << pathVids[i][j] << ", ";
    }

    std::cout << std::endl;

    for (int j = 0; j < pathFids.size(); j++) {
      std::cout << pathFids[j] << ", ";
    }

    std::cout << std::endl;
  }

  int originalVerticesNum = V.rows();
  V.conservativeResize(newVertices.size() + V.rows(), 3);
  for (int i = 0; i < newVertices.size(); i++) {
    for (int j = 0; j < 3; j++) {
      V(originalVerticesNum + i, j) = newVertices[i](j);
    }
  }
};
}
}
