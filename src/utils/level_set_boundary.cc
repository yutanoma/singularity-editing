#include "./level_set_boundary.h"

#include <array>
#include <iostream>

namespace rp {
namespace level_set_boundary {
void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &parameterization,
             const double &value, std::vector<Eigen::MatrixX3d> &paths) {
  Eigen::VectorXi hasBoundary(EF.rows());
  hasBoundary.setZero();

  for (int i = 0; i < EF.rows(); i++) {
    if (((parameterization(EV(i, 0)) - value) * (parameterization(EV(i, 1)) - value)) <= 0) {
      hasBoundary(i) = 1;
    }
  }

  Eigen::VectorXi isVisited(EF.rows());
  isVisited.setZero();

  while (hasBoundary.sum() > isVisited.sum()) {
    int smallestEdgeId = -1;
    for (int i = 0; i < hasBoundary.size(); i++) {
      if (hasBoundary(i) && !isVisited(i)) {
        smallestEdgeId = i;
        isVisited(i) = 1;
        break;
      }
    }

    if (smallestEdgeId == -1) {
      break;
    }

    int primalEdgeId = smallestEdgeId;

    std::array<std::vector<int>, 2> edgeIdPaths = {};

    for (int i = 0; i < 2; i++) {
      int currentFaceId = EF(primalEdgeId, i);
      int currentEdgeId = primalEdgeId;
      edgeIdPaths[i] = {currentEdgeId};

      if (currentFaceId == -1) {
        continue;
      }

      while (true) {
        int prevEdgeId = currentEdgeId;
        bool isEmptyEdge = false;

        for (int j = 0; j < 3; j++) {
          int eid = FE(currentFaceId, j);
          if (eid != prevEdgeId && hasBoundary(eid)) {
            currentEdgeId = eid;
            currentFaceId = EF(eid, 0) == currentFaceId ? EF(eid, 1) : EF(eid, 0);
            edgeIdPaths[i].emplace_back(eid);
            isVisited(eid) = 1;
            break;
          }
          if (j == 2) {
            isEmptyEdge = true;
          }
        }

        if (currentFaceId == -1) {
          break;
        } else {
          assert(!isEmptyEdge);
          if (currentEdgeId == primalEdgeId) {
            // 一周したということなのでループを終わらせる
            i += 1;
            break;
          }
        }
      }
    }

    std::vector<int> pathEdgeIds = {};

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < edgeIdPaths[i].size(); j++) {
        // iが0なら逆向きに挿入する
        int index = i == 0 ? edgeIdPaths[i].size() - 1 - j : j;
        if (i == 0 && edgeIdPaths[i + 1].size() > 0 && edgeIdPaths[i][edgeIdPaths[i].size() - 1] != edgeIdPaths[i][index]) {
          pathEdgeIds.emplace_back(edgeIdPaths[i][index]);
        }
      }
    }

    Eigen::MatrixX3d path(pathEdgeIds.size(), 3);
    for (int i = 0; i < pathEdgeIds.size(); i++) {
      int eid = pathEdgeIds[i];
      int baseVid = EV(eid, 0);
      int nextVid = EV(eid, 1);

      double baseVal = std::abs(parameterization(baseVid) - value);
      double nextVal = std::abs(parameterization(nextVid) - value);
      path.row(i) = (nextVal * V.row(baseVid) + baseVal * V.row(nextVid)) / (nextVal + baseVal);
    }

    paths.emplace_back(path);
  }

  std::cout << "paths num: " << paths.size() << std::endl;
}
}
}
