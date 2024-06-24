#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
namespace store {
namespace temporary_paths {
// outflow pathには流出するパスがすべて含まれる
void append_outflow_paths(const int &index, const Eigen::MatrixX3d &newPaths, const int &vid);

// inflow pathには流入するパスがすべて含まれる
void append_inflow_paths(const int &index, const Eigen::MatrixX3d &newPaths, const int &vid);

void reset_paths();

std::vector<Eigen::MatrixX3d> get_inflow_paths();
std::vector<Eigen::MatrixX3d> get_outflow_paths();

int get_last_inflow_vid();
int get_last_outflow_vid();
}
}
}
