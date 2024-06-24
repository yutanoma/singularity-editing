#pragma once

#include <Eigen/Core>
#include <vector>

namespace rp {
double parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                          const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                          const std::vector<int> &face_path,
                          const double &initial_angle);

double parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                          const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                          const std::vector<int> &face_path,
                          const double &initial_angle,
                          std::vector<int> &edgeIds,
                          std::vector<int> &edgeSigns);

double edges_parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                                const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                                const std::vector<int> &edge_path, const std::vector<int> &edge_signs,
                                const double &initial_angle);
}
