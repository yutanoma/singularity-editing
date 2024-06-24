#pragma once

#include <igl/remove_duplicate_vertices.h>

#include <Eigen/Core>

namespace rp {
void remove_duplicate_vertices(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F);

void remove_duplicate_vertices(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F, Eigen::VectorXi &J, const double &val);
}
