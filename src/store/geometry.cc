#include "./geometry.h"

#include <random>
#include <iostream>
#include <chrono>
#include <cstdlib>

#include <Eigen/Dense>

#include <igl/edge_topology.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/boundary_loop.h>
#include <igl/local_basis.h>
#include <igl/is_edge_manifold.h>
#include <igl/per_vertex_normals.h>
#include <directional/tree.h>

#include "../utils/vertex_edge_angles.h"
#include "../utils/hodge_star.h"
#include "../utils/open_path.h"
#include "../utils/non_dual_cycles.h"
#include "../utils/face_angles.h"
#include "../utils/curl_correction.h"
#include "../utils/stripe_patterns.h"
#include "../utils/index_prescription.h"
#include "../utils/rotation_to_parameterization.h"
#include "../utils/rotation_to_scalar.h"
#include "../utils/length_defects.h"
#include "../utils/angle_defects.h"
#include "../utils/angle_jumps.h"
#include "../utils/rotate_vector_field.h"
#include "../utils/parameterization_to_raw.h"
#include "../utils/singularity_decomposition.h"
#include "../utils/constrained_angles_cycle.h"
#include "../utils/constrained_parameters_cycle.h"
#include "../utils/convert_zero_vertex.h"
#include "../utils/fe_signs.h"
#include "../utils/angle_level_set_boundary.h"

namespace rp {
namespace store {

namespace {
Eigen::MatrixX3i initialF(0, 3);
Eigen::MatrixX3d initialV(0, 3);

double levelsetVal = .0;

// constant
Eigen::VectorXi innerEdges;

Eigen::MatrixX3i F(0, 3);
Eigen::MatrixX3d V(0, 3);
Eigen::MatrixXi EF(0, 2);
Eigen::MatrixXi FE(0, 2);
Eigen::MatrixXi EV(0, 2);
Eigen::MatrixXi FESigns(0, 2);

Eigen::SparseMatrix<double> identity(0, 0);

Eigen::MatrixX3d B1, B2, B3;

std::vector<std::vector<int>> VE = {};
std::vector<std::vector<double>> VEAngles = {};

Eigen::MatrixX3d faceAngles(0, 3);
Eigen::VectorXd anglesSum(0);

Eigen::VectorXd doubleAreas(0);

Eigen::MatrixXd normals(0, 3);

Eigen::VectorXi treeEdges(0), treeFathers(0);

Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> oneForm2ZeroFormSolver;
Eigen::SparseMatrix<double> zeroForm2TreeEdges;
bool needsFactorizationZeroform = true;

// triplets containing only basis loops and indexPresciription loops
std::vector<Eigen::Triplet<double>> basisLoops = {};
Eigen::SparseMatrix<double> basisLoopsMatrix;
int basisLoopsNum = 0;

// angle defect per loop
Eigen::VectorXd angleDefects(0);
Eigen::VectorXd angleJumps(0);

Eigen::SparseMatrix<double> hodgeStar;

std::vector<std::vector<int>> boundaryLoops;

// from here is related to vector field
// vector field assigned on cutV/cutF
Eigen::MatrixX3d vectorField(0, 3);
Eigen::VectorXd vectorAnglesPerVertex(0);

// parameterization defect
Eigen::VectorXd projectionPerEdge(0);

// rescalings
Eigen::VectorXd rescalings(0);

// the angle between the stripe pattern and the vectorfield. default is PI/2
double parameterizationAngle = igl::PI / 2;
Eigen::VectorXd rotatedVectorAnglesPerVertex(0);

// from here is related to the angled value stripe pattern
bool parameterizationComputationReady = false;
double velocity = 1;
double parameterGlobalRotation = 0;

// paths for the brushes
std::vector<std::vector<EdgePath>> brushPaths = {};

// the indices of the paths
std::vector<std::vector<EdgePath>> indexPrescriptionPaths = {};
std::vector<double> indexPrescriptionPathIndices = {};

Eigen::VectorXd parameterizationIndices(0);
Eigen::SparseMatrix<double> parameterizationMatrix;

Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> parameterizationSolver;

bool needsFactorizationParameterization = true;

// the 0-form of the stripe patterns
Eigen::VectorXd parameterization(0);
// 1-form
Eigen::VectorXd parameterizationOneForm(0);
// the levelset paths
std::vector<Eigen::MatrixX3d> paths = {};
};

void reset_parameterization() {
  parameterGlobalRotation = 0;

  parameterizationIndices.resize(basisLoopsNum);
  parameterizationIndices.setConstant(0);

  parameterizationMatrix = Eigen::SparseMatrix<double>();

  needsFactorizationParameterization = true;
}

void reset_geometry() {
  rp::fe_signs(F, FE, EV, FESigns);
  igl::local_basis(V, F, B1, B2, B3);
  rp::vertex_edge_angles::process(V, F, EV, EF, FE, VE, VEAngles, anglesSum);
  igl::doublearea(V, F, doubleAreas);

  boundaryLoops = {};
  igl::boundary_loop(F, boundaryLoops);

  igl::per_vertex_normals(V, F, normals);

  rp::non_dual_cycles::process(V, F, VE, FE, EF, EV, boundaryLoops, basisLoopsMatrix, basisLoops);
  basisLoopsNum = basisLoopsMatrix.rows();
  basisLoopsMatrix = Eigen::SparseMatrix<double>(basisLoopsNum, EV.rows());
  basisLoopsMatrix.setFromTriplets(basisLoops.begin(), basisLoops.end());
  rp::face_angles(V, F, B3, faceAngles);
  rp::angle_jumps(F, EV, VE, VEAngles, angleJumps);
  rp::angle_defects(F, V, anglesSum, faceAngles, basisLoopsMatrix, boundaryLoops, angleJumps, angleDefects);
  needsFactorizationZeroform = true;

  directional::tree(EV, treeEdges, treeFathers);

  identity.resize(EV.rows(), EV.rows());
  identity.setIdentity();

  innerEdges.resize(EV.rows());
  for (int i = 0; i < EV.rows(); i++)
    innerEdges(i) = i;

  rp::hodge_star(V, F, EF, EV, innerEdges, hodgeStar);

  rescalings.resize(V.rows());
  rescalings.setConstant(1);

  std::cout << "reset geometry finished" << std::endl;
}

void set_vectorfield(Eigen::MatrixXd &vf) {
  assert(vf.rows() == V.rows() && vf.cols() == 3);

  vectorField = vf;

  rescalings = Eigen::VectorXd::Ones(V.rows());

  vectorAnglesPerVertex.resize(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    int e = VE[i][0];
    int v0 = EV(e, 0), v1 = EV(e, 1);
    int otherV = v0 == i ? v1 : v0;

    Eigen::Vector3d n = normals.row(i).transpose().normalized();
    Eigen::Vector3d e0 = (V.row(otherV) - V.row(i)).transpose().normalized();
    e0 = e0 - n * n.dot(e0);

    auto e1 = n.cross(e0);

    double angle = std::atan2(e1.dot(vectorField.row(i)), e0.dot(vectorField.row(i)));

    vectorAnglesPerVertex(i) = angle;
  }

  rotate_vector_field(vectorAnglesPerVertex, parameterizationAngle, rotatedVectorAnglesPerVertex);

  rp::length_defects(V, F, EF, FE, EV, VE, VEAngles, rotatedVectorAnglesPerVertex, rescalings, velocity, projectionPerEdge);

  parameterizationComputationReady = true;
}

void reset_all() {
  V = initialV;
  F = initialF;

  igl::edge_topology(V, F, EV, FE, EF);

  reset_geometry();
}

void initialize(Eigen::MatrixX3i &_F, Eigen::MatrixX3d &_V) {
  initialF = _F;
  initialV = _V;

  V = _V;
  F = _F;

  igl::edge_topology(V, F, EV, FE, EF);

  reset_geometry();
};

void initialize(
  Eigen::MatrixX3i &_F,
  Eigen::MatrixX3d &_V,
  Eigen::MatrixXi &_EV,
  Eigen::MatrixXi &_FE,
  Eigen::MatrixXi &_EF
) {
  initialF = _F;
  initialV = _V;

  V = _V;
  F = _F;
  EV = _EV;
  FE = _FE;
  EF = _EF;

  reset_geometry();
}

void recompute_parameterization() {
  if (!parameterizationComputationReady) {
    return;
  }

  std::chrono::system_clock::time_point start, end;
  start = std::chrono::system_clock::now();
  bool factorized = false;

  if (needsFactorizationParameterization) {
    rp::constrained_parameters_cycle(V, F, EV, EF, treeFathers, projectionPerEdge, basisLoops, basisLoopsMatrix, indexPrescriptionPaths, indexPrescriptionPathIndices, brushPaths, 1, identity, levelsetVal, parameterGlobalRotation, parameterizationIndices, parameterizationMatrix, parameterizationSolver);

    needsFactorizationParameterization = false;
    factorized = true;
  }

  Eigen::VectorXd rotationAngle;
  double linfError;
  Eigen::SparseMatrix<double> identity(EV.rows(), EV.rows());
  identity.setIdentity();

  Eigen::VectorXd parameterizationDefects = parameterizationMatrix * projectionPerEdge;

  rp::index_prescription(V, F, EV, EF, innerEdges, parameterizationMatrix, parameterizationDefects, parameterizationIndices, parameterizationSolver, 1, identity, false, rotationAngle, linfError);
  parameterizationOneForm = rotationAngle + projectionPerEdge;
  rp::rotation_to_scalar(V, EV, parameterizationOneForm, .0, treeFathers, true, zeroForm2TreeEdges, oneForm2ZeroFormSolver, needsFactorizationZeroform, parameterization);
  // rp::rotation_to_parameterization(V, EV, parameterizationOneForm, .0, parameterization);

  end = std::chrono::system_clock::now();

  if (factorized) {
    std::cout << "[x] constrained_parameters_cycle: ";
  } else {
    std::cout << "[ ] constrained_parameters_cycle: ";
  }
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;
}

void level_set_boundary() {
  rp::angle_level_set_boundary::process(V, F, FE, EF, EV, parameterization, parameterizationOneForm, FESigns, parameterGlobalRotation,  paths);
}
  
void prescribe_parameterization_indices(const Eigen::VectorXd &indices) {
  Eigen::VectorXd _indices = indices;

  std::cout << indices.rows() << ", " << parameterizationIndices.rows() << ", " << basisLoopsNum + indexPrescriptionPaths.size() << std::endl;

  // only allow basis loops and indexPrescription loops to be specified
  // the order is basis loops first, then indexPrescription loops, then others
  if (_indices.rows() > basisLoopsNum + indexPrescriptionPaths.size()) {
    _indices = indices.head(basisLoopsNum + indexPrescriptionPaths.size());
  }

  std::cout << "l379" << std::endl;

  if (indices.head(F.rows() + boundaryLoops.size()).sum() != 0) {
    // must be zero
    std::cout << "warn: indices are invalid" << std::endl;
    return;
  }

  parameterizationIndices.head(indices.rows()) = indices;
  recompute_parameterization();
}

void set_parameterization_rotation(double &angle) {
  if (!parameterizationComputationReady) {
    std::cout << "Please set the vector field first" << std::endl;
    return;
  }

  parameterizationAngle = angle;
  rp::rotate_vector_field(vectorAnglesPerVertex, parameterizationAngle, rotatedVectorAnglesPerVertex);
  rp::length_defects(V, F, EF, FE, EV, VE, VEAngles, rotatedVectorAnglesPerVertex, rescalings, velocity, projectionPerEdge);
  recompute_parameterization();
}

void compute_rescalings() {
  if (!parameterizationComputationReady) {
    std::cout << "Please set the vector field first" << std::endl;
    return;
  }

  rp::curl_correction::process(V, F, FE, B1, B2, B3, VE, VEAngles, doubleAreas, EF, EV, 
                               faceAngles, vectorAnglesPerVertex, anglesSum, rescalings);
  recompute_parameterization();
}

void reset_rescalings() {
  if (!parameterizationComputationReady) {
    std::cout << "Please set the vector field first" << std::endl;
    return;
  }

  rescalings.setConstant(rescalings.rows(), 1);
  recompute_parameterization();
}

void set_default_stripe_patterns() {
  if (!parameterizationComputationReady) {
    std::cout << "Please set the vector field first" << std::endl;
    return;
  }

  reset_parameterization();

  rp::stripe_patterns::process(V, F, VE, VEAngles, doubleAreas, EF, EV, faceAngles, rotatedVectorAnglesPerVertex, rescalings, velocity, parameterization);
  rp::singularity_decomposition(EV, projectionPerEdge, parameterization, basisLoopsMatrix, parameterizationIndices);

  needsFactorizationParameterization = true;

  recompute_parameterization();
}

int get_brush_rows_num() {
  int num = 0;

  for (int i = 0; i < brushPaths.size(); i++) {
    num += brushPaths[i].size();
  }
  return num;
}

void set_brush_paths(std::vector<EdgePath> &brushPath) {
  if (!parameterizationComputationReady) {
    return;
  }

  brushPaths.emplace_back(brushPath);

  needsFactorizationParameterization = true;
}

void add_brush_path(const EdgePath &point, Eigen::MatrixX3d &newPath) {
  bool addNewPath = false;
  if (brushPaths.size() == 0 || brushPaths[brushPaths.size() - 1].size() == 0) {
    addNewPath = true;
  } else {
    auto lastpath = brushPaths[brushPaths.size() - 1];
    auto lastel = lastpath[lastpath.size() - 1];
    if (lastel.edgeId == point.edgeId) {
      newPath.resize(0, 3);
      return;
    }

    newPath.resize(2, 3);

    newPath.row(0) << V.row(EV(lastel.edgeId, 0)) + lastel.ratio * (V.row(EV(lastel.edgeId, 1)) - V.row(EV(lastel.edgeId, 0)));
    newPath.row(1) << V.row(EV(point.edgeId, 0)) + point.ratio * (V.row(EV(point.edgeId, 1)) - V.row(EV(point.edgeId, 0)));

    if (EV(lastel.edgeId, 0) == EV(point.edgeId, 0) || EV(lastel.edgeId, 0) == EV(point.edgeId, 1) || EV(lastel.edgeId, 1) == EV(point.edgeId, 0) || EV(lastel.edgeId, 1) == EV(point.edgeId, 1)) {
      addNewPath = false;
    } else {
      addNewPath = true;
    }
  }

  if (addNewPath) {
    std::vector<EdgePath> np = {point};
    brushPaths.emplace_back(np);
    newPath.resize(0, 3);
  } else {
    brushPaths[brushPaths.size() - 1].emplace_back(point);
  }

  needsFactorizationParameterization = true;
}

void set_index_prescription_paths(std::vector<EdgePath> &indexPrescriptionPath, int &index) {
  if (!parameterizationComputationReady) {
    return;
  }

  indexPrescriptionPaths.emplace_back(indexPrescriptionPath);
  indexPrescriptionPathIndices.emplace_back((double)index);

  needsFactorizationParameterization = true;
}

void set_prescription_path_index(int &pathId, int &index) {
  indexPrescriptionPathIndices[pathId] = (double)index;
}

void set_velocity(const double &v) {
  velocity = v;
  recompute_parameterization();
}

Eigen::MatrixX3i get_initial_f() {
  return initialF;
};
Eigen::MatrixX3d get_initial_v() {
  return initialV;
};
Eigen::MatrixX3i get_f() {
  return F;
};
Eigen::MatrixX3d get_v() {
  return V;
};
Eigen::MatrixXi get_ev() {
  return EV;
}
Eigen::MatrixXi get_fe() {
  return FE;
}
std::vector<std::vector<int>> get_ve() {
  return VE;
}
Eigen::VectorXd get_parameterization_indices() {
  return parameterizationIndices;
}
std::vector<std::vector<rp::EdgePath>> get_brush_paths() {
  return brushPaths;
}
std::vector<std::vector<rp::EdgePath>> get_index_prescription_paths() {
  return indexPrescriptionPaths;
}
std::vector<double> get_index_prescription_path_indices() {
  return indexPrescriptionPathIndices;
}
Eigen::VectorXd get_parameterization_diff() {
  return parameterizationOneForm;
}
Eigen::VectorXd get_parameterization() {
  return parameterization;
};
std::vector<Eigen::MatrixX3d> get_level_set_boundary() {
  return paths;
};
double get_velocity() {
  return velocity;
}
double get_rotation() {
  return parameterizationAngle;
}
Eigen::MatrixXd get_vector_field() {
  return vectorField;
};
};
}
