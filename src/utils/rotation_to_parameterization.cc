#include "rotation_to_parameterization.h"

#include <iostream>

namespace rp {
void rotation_to_parameterization(const Eigen::MatrixX3d &V, const Eigen::MatrixXi &EV, const Eigen::VectorXd &oneForm,
                                  const double &globalRotation,
                                  Eigen::VectorXd &periodicAngle) {
  // 1-formはEV(i, 1)-EV(i, 0)の値が入っている
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> ldltSolver;

  Eigen::SparseMatrix<std::complex<double>> d0(EV.rows(), V.rows() - 1), d0All(EV.rows(), V.rows());
  std::vector<Eigen::Triplet<std::complex<double>>> d0Triplets = {}, d0AllTriplets = {};
  for (int i = 0; i < EV.rows(); i++) {
    d0AllTriplets.emplace_back(Eigen::Triplet<std::complex<double>>(i, EV(i, 1), 1.));
    d0AllTriplets.emplace_back(Eigen::Triplet<std::complex<double>>(i, EV(i, 0), -std::exp(std::complex<double>(0, oneForm(i)))));

    if (EV(i, 1) != 0)
      d0Triplets.emplace_back(Eigen::Triplet<std::complex<double>>(i, EV(i, 1) - 1, 1.));
    if (EV(i, 0) != 0)
      d0Triplets.emplace_back(Eigen::Triplet<std::complex<double>>(i, EV(i, 0) - 1, -std::exp(std::complex<double>(0, oneForm(i)))));
  }

  d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
  d0All.setFromTriplets(d0AllTriplets.begin(), d0AllTriplets.end());

  Eigen::VectorXcd torhs = Eigen::VectorXcd::Zero(V.rows());
  torhs(0) = std::exp(std::complex<double>(0, globalRotation));  // global rotation
  Eigen::VectorXcd rhs = -d0All * torhs;
  ldltSolver.compute(d0.adjoint() * d0);

  assert(ldltSolver.info() == Eigen::Success);

  std::cout << "rotation_to_parmeterization solver computed" << std::endl;

  Eigen::VectorXcd periodicParams(V.rows());
  periodicParams(0) = std::exp(std::complex<double>(0, globalRotation));
  periodicParams.tail(V.rows() - 1) = ldltSolver.solve(d0.adjoint() * rhs);

  periodicAngle.resize(V.rows());

  for (int i = 0; i < periodicParams.rows(); i++) {
    periodicAngle(i) = std::atan2(periodicParams(i).imag(), periodicParams(i).real());
  }

  Eigen::VectorXcd diffVect = d0 * periodicParams.tail(V.rows() - 1) - rhs;

  std::cout << "rotation_to_parameterization linfError: " << diffVect.norm() << std::endl;

  // for (int i = 0; i < diffVect.rows(); i++) {
  //   double norm = std::sqrt(diffVect(i).real() * diffVect(i).real() + diffVect(i).imag() * diffVect(i).imag());
  //   if (norm > 0.001) {
  //     std::cout << "[" << i << "]: " << diffVect(i) << ", " << diffVect(i).real() * diffVect(i).real() + diffVect(i).imag() * diffVect(i).imag() << std::endl;
  //   }
  // }
}
}
