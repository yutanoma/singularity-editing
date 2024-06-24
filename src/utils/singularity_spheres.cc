#include "./singularity_spheres.h"

#include <directional/point_spheres.h>
#include <vector>
#include <iostream>

namespace rp {
  void singularity_spheres(const Eigen::MatrixX3i &F, const Eigen::MatrixX3d &V,
                           const Eigen::VectorXd &indices, Eigen::MatrixXi &sF, 
                           Eigen::MatrixXd &sV, Eigen::MatrixXd &sC, double radius) {
    std::vector<int> singularityFaceIds = {};
    std::vector<double> singularityIndices = {};
    int count = 0;

    for (int i = 0; i < F.rows(); i++) {
      if (indices(i) != 0) {
        singularityFaceIds.emplace_back(i);
        singularityIndices.emplace_back(indices(i));
        count++;
      }
    }

    std::cout << "singularity count: " << count << ", radius " << radius << std::endl;

    Eigen::MatrixX3d points(singularityFaceIds.size(), 3), C(singularityFaceIds.size(), 3);
    points.setZero();
    C.setZero();

    for (int i = 0; i < singularityFaceIds.size(); i++) {
      for (int j = 0; j < 3; j++) {
        points.row(i) += V.row(F(singularityFaceIds[i], j));
      }
      points.row(i) /= 3;

      if (singularityIndices[i] > 0) {
        C.row(i) << 1, 0, 0;
      } else {
        C.row(i) << 0, 0, 1;
      }
    }

    directional::point_spheres(points, radius, C, 10, sV, sF, sC);
  }
}