#include "./angle_level_set_boundary.h"

#include "vector_utils.h"
#include "remove_duplicate_vertices.h"

#include <array>
#include <iostream>

#include <igl/PI.h>
#include <igl/triangle/cdt.h>
#include <igl/is_edge_manifold.h>

namespace rp {
namespace angle_level_set_boundary {
double EPSILON = .0001;

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &zeroForm,
             const Eigen::VectorXd &oneForm, const Eigen::MatrixXi &FESigns,
             const double &value, std::vector<Eigen::MatrixX3d> &paths) {
  Eigen::MatrixX3d NV;
  Eigen::MatrixX3i NF;

  Eigen::MatrixX3d B1(F.rows(), 3);
  Eigen::MatrixX3d B2(F.rows(), 3);
  Eigen::MatrixX3d B3(F.rows(), 3);

  Eigen::VectorXd parameterization(V.rows());
  parameterization.setZero();

  std::vector<std::vector<int>> pathIds = {};

  process(V, F, FE, EF, EV, zeroForm, oneForm, FESigns, B1, B2, B3, value, false, paths, pathIds, NV, NF, parameterization);
}

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
             const Eigen::MatrixXi &EV, const Eigen::VectorXd &zeroForm,
             const Eigen::VectorXd &oneForm, const Eigen::MatrixXi &FESigns,
             const Eigen::MatrixX3d &B1, const Eigen::MatrixX3d &B2,
             const Eigen::MatrixX3d &B3,
             const double &value, const bool &updateGeometry,
             std::vector<Eigen::MatrixX3d> &paths,
             std::vector<std::vector<int>> &pathIds,
             Eigen::MatrixX3d &NV, Eigen::MatrixX3i &NF, Eigen::VectorXd &otherParameterization) {
  // -π～πに丸められた値
  double _value = rp::round_pi(value);

  if (otherParameterization.rows() != V.rows()) {
    int originalSize = otherParameterization.rows();
    otherParameterization.conservativeResize(V.rows());
    for (int i = 0; i < V.rows() - originalSize; i++) {
      otherParameterization(i + originalSize) = 0;
    }
  }

  paths = {};
  pathIds = {};

  std::vector<std::vector<int>> vertexIdsOnEdges(EF.rows(), std::vector<int>());

  int frontierVid = -1;
  std::vector<Eigen::Vector3d> newVertices = {};
  std::vector<int> newVerticesParams = {};

  for (int i = 0; i < EV.rows(); i++) {
    // 各edgeに頂点を記す
    int v0 = EV(i, 0), v1 = EV(i, 1);
    double val = rp::round_pi(zeroForm(v0)), diff = oneForm(i);
    double boundaryVal = value;
    int sign = diff > 0 ? 1 : -1;

    while ((boundaryVal - val) * sign < 0) {
      // boundaryValがvalより遅れていたら
      boundaryVal += sign * 2 * igl::PI;
    }

    std::vector<double> ratios = {};

    while (true) {
      double ratio = (boundaryVal - val) / diff;

      if (ratio > 1) {
        break;
      }

      if (ratio < 0) {
        boundaryVal += sign * 2 * igl::PI;
        continue;
      }

      if (ratio == 0) {
        ratio += .000001;
      } else if (ratio == 1) {
        ratio -= .000001;
      }

      ratios.emplace_back(ratio);
      boundaryVal += sign * 2 * igl::PI;
    }

    vertexIdsOnEdges[i].resize(ratios.size());
    for (int j = 0; j < ratios.size(); j++) {
      frontierVid += 1;
      vertexIdsOnEdges[i][j] = frontierVid;
      newVertices.emplace_back((V.row(v0) * (1 - ratios[j]) + V.row(v1) * ratios[j]).transpose());
      newVerticesParams.emplace_back(otherParameterization(v0) * (1 - ratios[j]) + otherParameterization(v1) * ratios[j]);
    }
  }

  std::vector<std::pair<int, int>> nextVertex(newVertices.size(), std::pair<int, int>(-1, -1));
  std::vector<std::vector<std::pair<int, int>>> edgePairsPerFace;

  if (updateGeometry) {
    edgePairsPerFace.resize(F.rows(), {});
  }

  for (int i = 0; i < F.rows(); i++) {
    double size = 0;
    std::array<std::vector<int>, 3> points({});
    for (int j = 0; j < 3; j++) {
      int edgeId = FE(i, j);
      points[j] = vertexIdsOnEdges[edgeId];

      // もしedgeとvertexが逆向きなら入れ替え
      int v0 = EV(edgeId, 0), v1 = EV(edgeId, 1);
      for (int k = 0; k < 3; k++) {
        if (F(i, k) != v0 && F(i, k) != v1) {
          if (F(i, (k + 1) % 3) != v0) {
            std::vector<int> vect = {};
            // std::cout << "l132" << std::endl;
            for (int l = 0; l < points[j].size(); l++) {
              // std::cout << l << ", " << points[j].size() << ", " << points[j].size() - l - 1 << ", " << std::endl;
              vect.emplace_back(points[j][points[j].size() - l - 1]);
            }
            points[j] = vect;
            // std::cout << "l136" << std::endl;
            break;
          }
        }
      }
    }

    std::array<int, 3> edgesOffset = {0, (int) points[0].size(), (int) points[0].size() + (int) points[1].size()};
    std::vector<int> connection(points[0].size() + points[1].size() + points[2].size(), -1);

    std::vector<int> realVid(connection.size(), -1);
    int count = 0;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < points[j].size(); k++) {
        realVid[count] = points[j][k];
        count++;
      }
    }

    // 両方向が下り坂or上り坂になっている場合に、頂点から近いものから順番につないでいく
    for (int j = 0; j < 3; j++) {
      int e0Id = FE(i, j);
      int e1Id = FE(i, (j + 1) % 3);

      double d0 = FESigns(i, j) * oneForm(e0Id);
      double d1 = FESigns(i, (j + 1) % 3) * oneForm(e1Id);

      if (d0 * d1 < 0) {
        int connectionCount = 0;
        // 下って上っている、もしくは上って下っている
        while (connectionCount / 2 < points[j].size() && connectionCount / 2 < points[(j + 1) % 3].size()) {
          int num = connectionCount / 2;
          int v0 = edgesOffset[j] + points[j].size() - num - 1;
          int v1 = edgesOffset[(j + 1) % 3] + num;
          connectionCount += 2;

          bool isV0Available = num < points[j].size();
          bool isV1Available = num < points[(j + 1) % 3].size();

          if (isV0Available && isV1Available) {
            if (connection[v0] == -1 && connection[v1] == -1) {
              connection[v0] = v1;
              connection[v1] = v0;
            } else if (connection[v0] != -1 || connection[v1] != -1) {
              // 関係者を全て洗いざらい-2にする
              for (int k = 0; k < connection.size(); k++) {
                if (connection[k] == v0 || connection[k] == v1) {
                  connection[k] = -2;
                }
              }

              connection[v0] = -2;
              connection[v1] = -2;
            }
          }
        }
      }
    }

    // -1か-2なら重心につなぐ
    int barycenterVid = frontierVid + 1;
    bool addBaryCenter = false;
    for (int j = 0; j < connection.size(); j++) {
      int vid = realVid[j];
      if (connection[j] == -1 || connection[j] == -2) {
        if (nextVertex[vid].first == -1) {
          nextVertex[vid].first = barycenterVid;
        } else {
          nextVertex[vid].second = barycenterVid;
        }
        addBaryCenter = true;
      } else {
        if (nextVertex[vid].first == -1) {
          nextVertex[vid].first = realVid[connection[j]];
        } else {
          nextVertex[vid].second = realVid[connection[j]];
        }
      }
    }

    if (addBaryCenter) {
      frontierVid += 1;
      newVertices.emplace_back(((V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3).transpose());
      newVerticesParams.emplace_back((otherParameterization(F(i, 0)) + otherParameterization(F(i, 1)) + otherParameterization(F(i, 2))) / 3);
      nextVertex.emplace_back(std::pair<int, int>(-1, -1));
    }

    // updateGeometryがtrueなら、分割するedgeをメモる
    if (updateGeometry && edgePairsPerFace.size() > i) {
      std::vector<bool> isConnected;
      isConnected.resize(connection.size(), false);
      for (int j = 0; j < connection.size(); j++) {
        if (!isConnected[j]) {
          isConnected[j] = true;

          int otherId = connection[j];

          if (otherId < 0) {
            edgePairsPerFace[i].emplace_back(std::pair<int, int>(realVid[j], barycenterVid));
          } else if (!isConnected[otherId]) {
            isConnected[otherId] = true;
            edgePairsPerFace[i].emplace_back(std::pair<int, int>(realVid[j], realVid[otherId]));
          }
        }
      }
    }
  }

  Eigen::VectorXi isVisited(newVertices.size());
  isVisited.setZero();

  while (newVertices.size() > isVisited.sum()) {
    int smallestVid = -1;
    for (int i = 0; i < newVertices.size(); i++) {
      if (!isVisited(i)) {
        smallestVid = i;
        isVisited(i) = 1;
        break;
      }
    }

    if (smallestVid == -1) {
      break;
    }

    std::array<std::vector<int>, 2> vidPaths = {};

    for (int i = 0; i < 2; i++) {
      vidPaths[i] = {};
      int currentVid = smallestVid;
      int nextVid = i == 0 ? nextVertex[currentVid].first : nextVertex[currentVid].second;

      vidPaths[i].emplace_back(currentVid);

      while (true) {
        if (nextVid == -1) {
          break;
        }

        vidPaths[i].emplace_back(nextVid);
        isVisited(nextVid) = 1;

        // ループをなしている場合
        if (nextVid == smallestVid) {
          i = 2;
          break;
        }

        int next = nextVertex[nextVid].first == currentVid
          ? nextVertex[nextVid].second
          : nextVertex[nextVid].first;

        currentVid = nextVid;
        nextVid = next;
      }
    }

    std::vector<int> pathVids = {};

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < vidPaths[i].size(); j++) {
        // iが0なら逆向きに挿入する
        int index = i == 0 ? vidPaths[i].size() - 1 - j : j;
        // iが1の時は、最初は飛ばす
        if (i == 0 || index != 0) {
          pathVids.emplace_back(vidPaths[i][index]);
        }
      }
    }

    Eigen::MatrixX3d path(pathVids.size(), 3);
    for (int i = 0; i < pathVids.size(); i++) {
      int vid = pathVids[i];
      path.row(i) << newVertices[vid].transpose();
    }

    paths.emplace_back(path);

    if (updateGeometry) {
      std::vector<int> realPathVids = {};
      realPathVids.resize(pathVids.size());
      for (int i = 0; i < pathVids.size(); i++) {
        auto vid = pathVids[i];
        realPathVids[i] = vid + V.rows();
      }

      // std::cout << "path[" << pathIds.size() << "]: ";
      // for (int i = 0; i < pathVids.size(); i++) {
      //   std::cout << realPathVids[i] << ", ";
      // }
      // std::cout << std::endl;
      pathIds.emplace_back(realPathVids);
    }
  }

  if (updateGeometry) {
    for (int i = 0; i < vertexIdsOnEdges.size(); i++) {
      for (int j = 0; j < vertexIdsOnEdges[i].size(); j++) {
        vertexIdsOnEdges[i][j] = vertexIdsOnEdges[i][j] + V.rows();
      }
    }

    int originalVerticesSize = V.rows();
    NV = V;
    NV.conservativeResize(originalVerticesSize + newVertices.size(), 3);
    otherParameterization.conservativeResize(originalVerticesSize + newVertices.size());
    for (int i = 0; i < newVertices.size(); i++) {
      NV.row(originalVerticesSize + i) = newVertices[i].transpose();
      otherParameterization(originalVerticesSize + i) = newVerticesParams[i];
    }

    Eigen::VectorXi isFaceUpdated(F.rows());
    isFaceUpdated.setZero();

    std::vector<std::array<int, 3>> newFaces = {};

    for (int i = 0; i < edgePairsPerFace.size(); i++) {
      if (edgePairsPerFace[i].size() > 0) {
        std::map<int, int> vid2localVid;
        std::vector<int> localVid2vid = {};

        Eigen::MatrixXi edges(edgePairsPerFace[i].size(), 2);
        for (int j = 0; j < edgePairsPerFace[i].size(); j++) {
          int v1 = edgePairsPerFace[i][j].first + originalVerticesSize;
          int v2 = edgePairsPerFace[i][j].second + originalVerticesSize;

          if (v1 > v2) {
            std::swap(v1, v2);
          }

          if (vid2localVid[v1] == 0 && (localVid2vid.size() == 0 || localVid2vid[0] != v1)) {
            vid2localVid[v1] = localVid2vid.size();
            localVid2vid.emplace_back(v1);
          }

          if (vid2localVid[v2] == 0 && (localVid2vid.size() == 0 || localVid2vid[0] != v2)) {
            vid2localVid[v2] = localVid2vid.size();
            localVid2vid.emplace_back(v2);
          }

          edges.row(j) << std::min(vid2localVid[v1], vid2localVid[v2]), std::max(vid2localVid[v1], vid2localVid[v2]);
        }

        bool debug = false;

        for (int j = 0; j < 3; j++) {
          int vid = F(i, j);
          vid2localVid[vid] = localVid2vid.size();
          vid2localVid[vid] = localVid2vid.size();
          localVid2vid.emplace_back(vid);

          // if (vid == 1362 || vid == 1365 || vid == 1366 || vid == 1343 || vid == 1371 || vid == 1375 || vid == 1375 || vid == 1381 || vid == 1387 || vid == 1389 || vid == 1403 ||  vid == 1427 || i == 3057 || i == 4021 || i == 4022 || i == 5536 || i == 5537 || i == 6415 || i == 6416 || i == 7155 || i == 7156 || i == 9198 || i == 11724 || i == 11725 || i == 12880 || i == 12881 || i == 13228 || i == 13472 || i == 14045 || i == 14420) {
          //   debug = true;
          // }
        }

        std::vector<int> vidsClockwise = {};

        for (int j = 0; j < 3; j++) {
          int edgeId = FE(i, j), sign = FESigns(i, j);
          // std::cout << edgeId << ", " << vertexIdsOnEdges.size() << std::endl;
          for (int k = 0; k < vertexIdsOnEdges[edgeId].size(); k++) {
            int idx = sign > 0 ? k : vertexIdsOnEdges[edgeId].size() - 1 - k;
            vidsClockwise.emplace_back(vertexIdsOnEdges[edgeId][idx]);
          }
          int endVid = sign > 0 ? EV(edgeId, 1) : EV(edgeId, 0);
          vidsClockwise.emplace_back(endVid);
        }

        edges.conservativeResize(edges.rows() + vidsClockwise.size(), 2);

        for (int j = 0; j < vidsClockwise.size(); j++) {
          int v0 = vidsClockwise[j];
          int v1 = vidsClockwise[(j + 1) % vidsClockwise.size()];
          edges(edgePairsPerFace[i].size() + j, 0) = vid2localVid[v0];
          edges(edgePairsPerFace[i].size() + j, 1) = vid2localVid[v1];
        }

        // if (edges.rows() > 2) {
        //   std::cout << std::endl;
        // }

        Eigen::MatrixXd vertices(localVid2vid.size(), 2);

        // std::cout << "face[" << i << "]: " << F.row(i) << std::endl;

        auto b1 = B1.row(i).transpose().normalized();
        auto b2 = B2.row(i).transpose().normalized();

        for (int j = 0; j < localVid2vid.size(); j++) {
          auto p = NV.row(localVid2vid[j]).transpose();
          vertices(j, 0) = b1.dot(p);
          vertices(j, 1) = b2.dot(p);
        }

        // std::cout << "edges: " << edges << std::endl;
        // std::cout << "vertices: " << vertices << std::endl;
        // verticesにはB1とB2への射影を表示
        Eigen::MatrixXd newV;
        Eigen::MatrixXi newF, newE;
        Eigen::VectorXi J;
        // igl::triangle::cdt(vertices, edges, "-c", newV, newF, newE, J);
        // std::cout << vertices << std::endl;
        // std::cout << std::endl;
        // std::cout << edges << std::endl;
        igl::triangle::triangulate(vertices, edges, Eigen::MatrixXi(0, 2), "-Q", newV, newF);

        if (debug) {
          std::cout << F.row(i) << std::endl;
        }

        if (newV.rows() != vertices.rows()) {
          std::cout << vertices << std::endl;
          std::cout << edges << std::endl;
          std::cout << newV << std::endl;
          std::cout << newF << std::endl;
        }

        assert(newV.rows() == vertices.rows());
        // assert(igl::is_edge_manifold(newF));

        // std::cout << newF << std::endl;

        Eigen::Vector3d originalNormal = B3.row(i).transpose().normalized();

        for (int j = 0; j < newF.rows(); j++) {
          int v0 = localVid2vid[newF(j, 0)];
          int v1 = localVid2vid[newF(j, 1)];
          int v2 = localVid2vid[newF(j, 2)];

          Eigen::Vector3d e0 = NV.row(v1).transpose() - NV.row(v0).transpose();
          Eigen::Vector3d e1 = NV.row(v2).transpose() - NV.row(v0).transpose();
          Eigen::Vector3d nn = e0.cross(e1).normalized();

          bool isInverted = nn.dot(originalNormal) < 0;

          bool localDebug = false;

          // 3点が同一直線上にあるかどうかを判定
          // bool isAvailable = true;
          // std::array<int, 3> vs = {v0, v1, v2};
          // for (int k = 0; k < 3; k++) {
          //   double d1 = (NV.row(vs[k]) - NV.row(vs[(k + 1) % 3])).norm();
          //   double d2 = (NV.row(vs[k]) - NV.row(vs[(k + 2) % 3])).norm();
          //   double d3 = (NV.row(vs[(k + 1) % 3]) - NV.row(vs[(k + 2) % 3])).norm();
          //   if ((d1 + d2 - d3) < .0001) {
          //     isAvailable = false;
          //   }
          // }

          // if (!isAvailable) {
          //   continue;
          // }

          std::array<int, 3> nf;
          for (int k = 0; k < 3; k++) {
            int vid = localVid2vid[newF(j, isInverted ? 2 - k : k)];
            nf[k] = vid;
            // if (vid == 1362 || vid == 1365 || vid == 1366 || vid == 1343 || vid == 1371 || vid == 1375 || vid == 1375 || vid == 1381 || vid == 1387 || vid == 1389 || vid == 1403 ||  vid == 1427) {
            //   localDebug = true;
            // }
          }

          if (debug || localDebug) {
            std::cout << "(" << nf[0] << " " << nf[1] << " " << nf[2] << "), inverted: " << isInverted << std::endl;
            std::cout << nn.transpose() << ", " << originalNormal.transpose() << ", " << nn.dot(originalNormal) << std::endl << std::endl;
          }

          newFaces.emplace_back(nf);
        }

        isFaceUpdated(i) = 1;
      }
    }

    NF.resize(F.rows() + newFaces.size() - isFaceUpdated.sum(), 3);

    int front = 0;
    for (int i = 0; i < F.rows(); i++) {
      if (!isFaceUpdated(i)) {
        NF.row(front) = F.row(i);
        front++;
      }
    }
    for (int i = 0; i < newFaces.size(); i++) {
      NF.row(front) << newFaces[i][0], newFaces[i][1], newFaces[i][2];
      front++;
    }

    // 最後、重複してるvertexを消す
    // Eigen::VectorXi J;
    // int os = NV.rows();
    // auto of = NF.rows();

    // std::cout << "505" << std::endl;
    // rp::remove_duplicate_vertices(NV, NF, J, EPSILON);
    // std::cout << "507" << std::endl;
    // std::cout << J.rows() << ", " << NV.rows() << ", " << os << ", " << NF.rows() << ", " << of << std::endl;
    // std::vector<std::vector<int>> _pathIds = {};
    // _pathIds.resize(pathIds.size(), std::vector<int>());
    // for (int i = 0; i < pathIds.size(); i++) {
    //   for (int j = 0; j < pathIds[i].size(); j++) {
    //     int id = J(pathIds[i][j]);
    //     if (_pathIds[i].size() == 0 || id != _pathIds[i][_pathIds[i].size() - 1]) {
    //       _pathIds[i].emplace_back(id);
    //     }
    //   }
    // }
    // pathIds = _pathIds;
    // std::cout << "521" << std::endl;
  }
}
}
}
