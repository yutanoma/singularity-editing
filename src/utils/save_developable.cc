#include "save_developable.h"

#include "divide_chunks.h"
#include "triangle_angles.h"
#include <igl/edge_topology.h>
#include <stack>
#include <map>
#include <array>
#include <vector>
#include <deque>
#include <fstream>
#include <iostream>

namespace rp {
  void save_developable(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const std::string &filename) {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F, EV, FE, EF);

    std::map<int, bool> unVisitedFaces;

    Eigen::VectorXi isVisited(F.rows());
    isVisited.setZero();

    int count = 0;

    while (isVisited.sum() < F.rows()) {
      int fid = -1;
      for (int i = 0; i < isVisited.rows(); i++) {
        if (isVisited(i) == 0) {
          fid = i;
        }
      }
      if (fid == -1) {
        break;
      }

      std::vector<int> nextFaces = {};

      for (int j = 0; j < 3; j++) {
        int edgeId = FE(fid, j);
        if (edgeId != -1) {
          int next = EF(edgeId, 0) == fid ? EF(edgeId, 1) : EF(edgeId, 0);
          if (next == -1) {
            continue;
          }
          if (isVisited(next) == 0) {
            nextFaces.emplace_back(next);
          }
        }
      }

      std::cout << "l63" << std::endl;
      isVisited(fid) = 1;

      std::deque<int> fidsList = {fid};
      int size = nextFaces.size();

      for (int i = 0; i < std::min(size, 2); i++) {
        int currentFid = nextFaces[i];
        int prevFid = fid;

        while (true) {
          if (i == 0) {
            fidsList.emplace_back(currentFid);
          } else {
            fidsList.emplace_front(currentFid);
          }

          // std::cout << currentFid << std::endl;

          isVisited(currentFid) = 1;

          std::vector<int> nf = {};

          for (int j = 0; j < 3; j++)
          {
            int edgeId = FE(currentFid, j);
            if (edgeId != -1)
            {
              if (EF(edgeId, 0) != prevFid && EF(edgeId, 1) != prevFid)
              {
                int next = EF(edgeId, 0) == currentFid ? EF(edgeId, 1) : EF(edgeId, 0);
                if (next == -1) {
                  continue;
                }
                if (isVisited(next) == 0)
                {
                  nf.emplace_back(next);
                }
              }
            }
          }

          if (nf.size() == 0) {
            break;
          }

          currentFid = nf[0];
          prevFid = currentFid;
        }
      }

      std::cout << "l101: " << fidsList.size() << std::endl;

      std::vector<std::array<std::array<double, 2>, 3>> coordList = {};
      coordList.resize(fidsList.size());

      for (int i = 0; i < fidsList.size(); i++) {
        if (i == 0) {
          double length = (V.row(F(fidsList[i], 1)) - V.row(F(fidsList[i], 0))).norm();
          double a = length;
          double b = (V.row(F(fidsList[i], 2)) - V.row(F(fidsList[i], 1))).norm();
          double c = (V.row(F(fidsList[i], 0)) - V.row(F(fidsList[i], 2))).norm();
          double A, B, C;
          rp::triangle_angles(a, b, c, A, B, C);

          std::array<std::array<double, 2>, 3> coord;
          coord[0] = {.0, .0};
          coord[1] = {.0, length};
          coord[2] = {c * std::cos(B), c * std::sin(B)};
          coordList[i] = coord;
        } else {
          // その三角形を前の三角形に接続する
          int prevFid = fidsList[i - 1];
          int currFid = fidsList[i];

          // 一致している辺を見つける
          int v0idx = -1, v1idx = -1, v0 = -1, v1 = -1, baseIdx = -1;
          for (int j = 0; j < 3; j++) {
            baseIdx = j;
            v0 = F(currFid, j), v1 = F(currFid, (j + 1) % 3);
            for (int k = 0; k < 3; k++) {
              if (F(prevFid, k) != v0 && F(prevFid, k) != v1 && F(prevFid, (k + 1) % 3) == v1 && F(prevFid, (k + 2) % 3) == v0) {
                // その次とそのあとが共有する辺になっている
                v1idx = (k + 1) % 3;
                v0idx = (k + 2) % 3;
 
                break;
              }
            }

            if (v0idx != -1) {
              break;
            }
          }

          assert(v0idx != -1 && v1idx != -1);

          std::array<double, 2> v0coord = coordList[i-1][v0idx];
          std::array<double, 2> v1coord = coordList[i-1][v1idx];

          double a = (V.row(F(currFid, (baseIdx + 1) % 3)) - V.row(F(currFid, baseIdx))).norm();
          double b = (V.row(F(currFid, (baseIdx + 2) % 3)) - V.row(F(currFid, (baseIdx + 1) % 3))).norm();
          double c = (V.row(F(currFid, baseIdx)) - V.row(F(currFid, (baseIdx + 2) % 3))).norm();
          double A, B, C;
          rp::triangle_angles(a, b, c, A, B, C);

          std::array<std::array<double, 2>, 3> coord;
          coord[0] = v0coord;
          coord[1] = v1coord;
          double x = v1coord[0] - v0coord[0];
          double y = v1coord[1] - v0coord[1];
          double nx = x * std::cos(B) - y * std::sin(B);
          double ny = x * std::sin(B) + y * std::cos(B);
          double norm = std::sqrt(nx * nx + ny * ny);
          nx = nx * c / norm;
          ny = ny * c / norm;

          coord[2] = {nx + v0coord[0], ny + v0coord[1]};
          coordList[i] = coord;
        }
      }

      std::ofstream writingFile;
      writingFile.open(filename + "_" + std::to_string(count) + ".svg", std::ios::trunc);
      count++;
      writingFile.clear();
      writingFile << "<?xml version=\"1.0\"?>" << std::endl;
      writingFile << "<svg xmlns=\"http://www.w3.org/2000/svg\">";
      for (int i = 0; i < coordList.size(); i++) {
        writingFile << "<polygon fill=\"black\" points=\"";
        for (int j = 0; j < coordList[i].size(); j++) {
          for (int k = 0; k < coordList[i][j].size(); k++) {
            writingFile << coordList[i][j][k];
            if (k == 0) {
              writingFile << ",";
            } else {
              writingFile << " ";
            }
          }
        }
        writingFile << "\" />";
      }
      writingFile << "</svg>";
      writingFile.close();
    }
  }
}