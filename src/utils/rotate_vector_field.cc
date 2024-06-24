#include "./rotate_vector_field.h"

namespace rp {
void rotate_vector_field(const Eigen::VectorXd &originalVectorAngles, const double &angle,
                          Eigen::VectorXd &rotatedVectorAngles) {
  rotatedVectorAngles.resize(originalVectorAngles.rows());

  for (int i = 0; i < originalVectorAngles.rows(); i++) {
    rotatedVectorAngles(i) = rp::round_pi(originalVectorAngles(i) + angle);
  }
};
}
