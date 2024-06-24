#include "./temporary_paths.h"

namespace rp {
namespace store {
namespace temporary_paths {
namespace {
  std::vector<Eigen::MatrixX3d> outflowPaths = {};
  std::vector<std::vector<int>> outflowPathVids = {};
  std::vector<Eigen::MatrixX3d> inflowPaths = {};
  std::vector<std::vector<int>> inflowPathVids = {};
}

// outflow pathには流出するパスがすべて含まれる
void append_outflow_paths(const int &index, const Eigen::MatrixX3d &newPaths, const int &vid) {
  if (outflowPaths.size() <= index) {
    outflowPaths.emplace_back(newPaths);
    std::vector<int> newVids = {vid};
    outflowPathVids.emplace_back(newVids);
  } else {
    int rowSize = outflowPaths[index].rows();
    outflowPaths[index].conservativeResize(rowSize + newPaths.rows() - 1, 3);

    for (int i = 0; i < newPaths.rows() - 1; i++) {
      for (int j = 0; j < 3; j++) {
        outflowPaths[index](rowSize + i, j) = newPaths(i + 1, j);
      }
    }

    outflowPathVids[index].emplace_back(vid);
  }
};

// inflow pathには流入するパスがすべて含まれる
void append_inflow_paths(const int &index, const Eigen::MatrixX3d &newPaths, const int &vid) {
  if (inflowPaths.size() <= index) {
    inflowPaths.emplace_back(newPaths);
    std::vector<int> newVids = {vid};
    inflowPathVids.emplace_back(newVids);
  } else {
    int rowSize = inflowPaths[index].rows();
    inflowPaths[index].conservativeResize(rowSize + newPaths.rows() - 1, 3);

    for (int i = 0; i < newPaths.rows() - 1; i++) {
      for (int j = 0; j < 3; j++) {
        inflowPaths[index](rowSize + i, j) = newPaths(i + 1, j);
      }
    }

    inflowPathVids[index].emplace_back(vid);
  }
}

void reset_paths() {
  inflowPaths = {};
  outflowPaths = {};
}

std::vector<Eigen::MatrixX3d> get_inflow_paths() {
  return inflowPaths;
};

std::vector<Eigen::MatrixX3d> get_outflow_paths() {
  return outflowPaths;
};

int get_last_inflow_vid() {
  if (inflowPathVids.size() == 0) {
    return -1;
  }
  if (inflowPathVids[inflowPathVids.size() - 1].size() == 0) {
    return -1;
  }

  auto lastPath = inflowPathVids[inflowPathVids.size() - 1];
  return lastPath[lastPath.size() - 1];
}

int get_last_outflow_vid() {
  if (outflowPathVids.size() == 0) {
    return -1;
  }
  if (outflowPathVids[outflowPathVids.size() - 1].size() == 0) {
    return -1;
  }

  auto lastPath = outflowPathVids[outflowPathVids.size() - 1];
  return lastPath[lastPath.size() - 1];
}
}
}
}
