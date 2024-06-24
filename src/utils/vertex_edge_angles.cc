#include "./vertex_edge_angles.h"

#include "./vector_utils.h"

#include <deque>
#include <iostream>
#include <igl/PI.h>

namespace rp {
namespace vertex_edge_angles {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi EV,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
             std::vector<std::vector<int>> &VE, std::vector<std::vector<double>> &VEAngles) {
  Eigen::VectorXd anglesSumPerVertex;
  process(V, F, EV, EF, FE, VE, VEAngles, anglesSumPerVertex);
}

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::MatrixXi EV,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
             std::vector<std::vector<int>> &VE, std::vector<std::vector<double>> &VEAngles, Eigen::VectorXd &anglesSumPerVertex) {
  VE.resize(V.rows(), {});
  VEAngles.resize(V.rows(), {});
  anglesSumPerVertex.resize(V.rows());
  anglesSumPerVertex.setConstant(2 * igl::PI);

  // まず仮のVEを作る。これはVEの中に順序がついていない。
  std::vector<std::vector<int>> temporaryVE;
  temporaryVE.resize(V.rows(), {});

  for (int i = 0; i < EV.rows(); i++) {
    for (int j = 0; j < EV.cols(); j++) {
      temporaryVE[EV(i, j)].emplace_back(i);
    }
  }

  // 次に、VEを作っていく。このVEは0番目の要素からedgeIdが反時計回りに順番に入っている
  for (int i = 0; i < temporaryVE.size(); i++) {
    if (temporaryVE[i].size() == 0) {
      continue;
    }

    int currentEdgeId = temporaryVE[i][0];
    std::deque<int> edgeIds = {currentEdgeId};
    // 1なら反時計回りに、-1なら時計回りにedgeを回収していく
    int iterationDirection = 1;

    while (edgeIds.size() < temporaryVE[i].size()) {
      auto fids = EF.row(currentEdgeId);
      int nextEdgeId = -1;
      for (int j = 0; j < 2; j++) {
        if (fids(j) == -1) {
          continue;
        }
        int fid = fids(j);
        for (int k = 0; k < 3; k++) {
          int eid = FE(fid, k);
          int nextEid = FE(fid, (k + 3 - iterationDirection) % 3);
          if (eid == currentEdgeId && (EV(nextEid, 0) == i || EV(nextEid, 1) == i)) {
            nextEdgeId = nextEid;
          }
        }
      }
      if (nextEdgeId == -1) {
        assert(iterationDirection == 1);
        iterationDirection = -iterationDirection;
        currentEdgeId = temporaryVE[i][0];
      } else {
        currentEdgeId = nextEdgeId;
        if (iterationDirection == 1) {
          edgeIds.emplace_back(nextEdgeId);
        } else if (iterationDirection == -1) {
          edgeIds.emplace_front(nextEdgeId);
        }
      }
    }

    VE[i].resize(edgeIds.size());
    for (int j = 0; j < edgeIds.size(); j++) {
      VE[i][j] = edgeIds[j];
    }
  }

  // 最後に、最初のangleから目的のedgeまでの角度を格納していく
  for (int i = 0; i < VE.size(); i++) {
    assert(VE[i].size() != 0);

    VEAngles[i].resize(VE[i].size(), 0);

    int firstEdge = VE[i][0], lastEdge = VE[i][VE[i].size() - 1];

    bool isHalfdisk = false;
    if (EF(firstEdge, 0) == -1 || EF(firstEdge, 1) == -1 || EF(lastEdge, 0) == -1 || EF(lastEdge, 1) == -1) {
      // このvertexの周りは閉じていない、半ディスク同相になっている
      isHalfdisk = true;
    }

    // edge i->edge i+1の角度をiに保管している
    std::vector<double> angles;
    angles.resize(VE[i].size());

    for (int j = 0; j < VE[i].size(); j++) {
      int eid = VE[i][j];
      int vid = EV(eid, 0) == i ? EV(eid, 1) : EV(eid, 0);
      assert(EV(eid, 0) == i || EV(eid, 1) == i);
      Eigen::Vector3d edge = (V.row(vid) - V.row(i)).transpose().normalized();

      int nextEid = VE[i][(j + 1) % VE[i].size()];
      int nextVid = EV(nextEid, 0) == i ? EV(nextEid, 1) : EV(nextEid, 0);
      Eigen::Vector3d nextEdge = (V.row(nextVid) - V.row(i)).transpose().normalized();
      assert(EV(nextEid, 0) == i || EV(nextEid, 1) == i);

      double dot = edge.dot(nextEdge);
      double angle = std::acos(dot);

      angles[j] = angle;
    }

    // halfdiskならangleの合計値はそのまま、diskならangleの合計値は2πに丸める
    double angleSum = .0;
    if (isHalfdisk) {
      angleSum = 2 * igl::PI;
    } else {
      for (auto a : angles) {
        angleSum += a;
      }
    }
    anglesSumPerVertex(i) = angleSum;

    VEAngles[i][0] = 0;
    for (int j = 0; j < VE[i].size() - 1; j++) {
      VEAngles[i][j + 1] = VEAngles[i][j] + angles[j] * 2 * igl::PI / angleSum;
    }
  }
}
}
}
