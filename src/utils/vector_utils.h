#pragma once

#include <igl/PI.h>
#include <vector>

#include <Eigen/Core>

namespace rp {
double get_angle(const Eigen::Vector3d &basis_vec,
                 const Eigen::Vector3d &angled_vec,
                 const Eigen::Vector3d &normal);

Eigen::Vector3d face_normal(const Eigen::MatrixX3d &V,
                            const Eigen::MatrixX3i &F, const int face_id);

void vector_decomposition(Eigen::Vector2d &x, Eigen::Vector2d &v_a,
                          Eigen::Vector2d &v_b,
                          Eigen::Vector2d &result);

void vector_decomposition_3d(Eigen::Vector3d &x, Eigen::Vector3d &v_a,
                             Eigen::Vector3d &v_b,
                             Eigen::Vector2d &result);

void get_edge_from_faces(const int &prevFid, const int &postFid,
                         const Eigen::MatrixXi &FE, const Eigen::MatrixXi &EF,
                         int &edge, int &edgeSign);

double round_pi(const double &angle);

void create_tree_from_tEf(
    const Eigen::VectorXi& dualTreeFathers, const Eigen::MatrixXi& EV,
    std::vector<std::vector<int>>& treeNextNodeList);
}  // namespace rp
