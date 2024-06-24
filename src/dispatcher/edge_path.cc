#include "./edge_path.h"

#include "../utils/exact_geodesic.h"
#include "../store/geometry.h"
#include "../store/temporary_paths.h"

#include <iostream>

namespace rp {
namespace dispatcher {
  void get_edge_path(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                     const Eigen::MatrixXi &EV, const std::vector<std::vector<int>> &VE,
                     const int &nextVid, const int &prevVid, Eigen::MatrixX3d &newPaths, std::vector<std::vector<EdgePath>> &paths) {
    Eigen::VectorXi positions(2);
    positions(0) = prevVid;
    positions(1) = nextVid;

    std::vector<Eigen::MatrixX3d> resultPoints;
    std::vector<Eigen::VectorXi> p = {positions};

    auto _V = V;
    auto _F = F;
    std::vector<Eigen::VectorXi> arr;
    std::vector<std::vector<std::pair<std::array<int, 2>, double>>> edgePaths;

    rp::exact_geodesic::split_path(_V, _F, p, false, resultPoints, arr, edgePaths);

    if (resultPoints.size() > 0) {
      newPaths = resultPoints[0];
    } else {
      newPaths.resize(0, 3);
    }

    // for (auto ep : edgePaths) {
    //   for (auto p : ep) {
    //     std::cout << "(" << p.first[0] << ", " << p.first[1] << ", " << p.second << "), ";
    //   }
    //   std::cout << std::endl;
    // }

    paths.resize(edgePaths.size(), {});
    for (int i = 0; i < edgePaths.size(); i++) {
      for (int j = 0; j < edgePaths[i].size(); j++) {
        auto original = edgePaths[i][j];
        int v1 = original.first[0];
        int v2 = original.first[1];
        double ratio = original.second;

        int edgeId = -1;
        for (int k = 0; k < VE[v1].size(); k++) {
          int eid = VE[v1][k];
          if (EV(eid, 0) == v2) {
            edgeId = eid;
            ratio = 1.0 - ratio;
          } else if (EV(eid, 1) == v2) {
            edgeId = eid;
          }
        }

        assert(edgeId != -1);

        paths[i].emplace_back(EdgePath(edgeId, ratio));
      }
    }
  }

  void set_brush_paths(const int &prevVid, const int &nextVid, Eigen::MatrixX3d &newPaths) {
    auto F = rp::store::get_f();
    auto V = rp::store::get_v();
    auto EV = rp::store::get_ev();
    auto VE = rp::store::get_ve();
    std::vector<std::vector<EdgePath>> edgePaths = {};

    get_edge_path(V, F, EV, VE, nextVid, prevVid, newPaths, edgePaths);

    for (auto ep : edgePaths) {
      rp::store::set_brush_paths(ep);
    }
  }

  void set_brush_stroke(const int &_v0, const int &_v1, const double &_ratio, Eigen::MatrixX3d &newPaths) {
    auto EV = rp::store::get_ev();
    auto VE = rp::store::get_ve();

    int v0 = std::min(_v0, _v1), v1 = std::max(_v0, _v1);
    double ratio = v0 == _v0 ? _ratio : 1 - _ratio;

    int edgeId = -1;
    std::cout << "l80: ";
    for (int i = 0; i < VE[v0].size(); i++) {
      std::cout << EV(VE[v0][i], 1) << ", ";
      if (EV(VE[v0][i], 1) == v1) {
        edgeId = VE[v0][i];
        break;
      }
    }
    std::cout << std::endl;
    if (edgeId == -1) {
      newPaths.resize(0, 3);
      return;
    }

    std::cout << "l91" << std::endl;

    rp::store::add_brush_path(EdgePath(edgeId, ratio), newPaths);

    std::cout << "l95" << std::endl;
  }

  void set_index_prescription_paths() {
    // TODO
  }

  void set_prescription_path_index() {
    // TODO
  }
};
};
