#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>

#include "../utils/edge_paths.h"
#include "../utils/constrained_parameters_cycle.h"

namespace rp {
namespace store {
void initialize(Eigen::MatrixX3i &_F, Eigen::MatrixX3d &_V);
void initialize(Eigen::MatrixX3i &_F, Eigen::MatrixX3d &_V, Eigen::MatrixXi &_EV, Eigen::MatrixXi &_FE, Eigen::MatrixXi &_EF);

void set_vectorfield(Eigen::MatrixXd &vf);

void reset_parameterization();
void reset_all();
void recompute_parameterization();
void prescribe_parameterization_indices(const Eigen::VectorXd &indices);
void set_parameterization_rotation(double &angle);
void compute_rescalings();
void level_set_boundary();
void reset_rescalings();
void set_default_stripe_patterns();
void set_brush_paths(std::vector<EdgePath> &brushPath);
void add_brush_path(const EdgePath &point, Eigen::MatrixX3d &newPath);
void set_index_prescription_paths(std::vector<EdgePath> &indexPrescriptionPath, int &index);
void set_prescription_path_index(int &pathId, int &index);
void set_velocity(const double &v);

Eigen::MatrixX3i get_initial_f();
Eigen::MatrixX3d get_initial_v();
Eigen::MatrixX3i get_f();
Eigen::MatrixX3d get_v();
Eigen::MatrixXi get_ev();
Eigen::MatrixXi get_fe();
std::vector<std::vector<int>> get_ve();
Eigen::VectorXd get_parameterization_indices();
std::vector<std::vector<rp::EdgePath>> get_brush_paths();
std::vector<std::vector<rp::EdgePath>> get_index_prescription_paths();
std::vector<double> get_index_prescription_path_indices();
Eigen::VectorXd get_parameterization();
Eigen::VectorXd get_parameterization_diff();
double get_rotation();
double get_velocity();
std::vector<Eigen::MatrixX3d> get_level_set_boundary();
Eigen::MatrixXd get_vector_field();
Eigen::MatrixXd get_perpendicular_vector_field();
};
};
