#include "angle_color.h"

#include <igl/PI.h>
#include "vector_utils.h"

namespace rp {
  void get_rgb_from_angle(const double &angle, double &r, double &g, double &b) {
    if (0 <= angle && angle < igl::PI / 3) {
      r = 1;
      g = 1. * angle / (igl::PI / 3);
      b = 0;
    } else if (igl::PI / 3 <= angle && angle < igl::PI * 2 / 3) {
      r = 1. * (igl::PI * 2 / 3 - angle) / (igl::PI / 3);
      g = 1;
      b = 0;
    } else if (igl::PI * 2 / 3 <= angle && angle < igl::PI) {
      r = 0;
      g = 1;
      b = 1. * (angle - igl::PI * 2 / 3) / (igl::PI / 3);
    } else if (igl::PI <= angle && angle < igl::PI * 4 / 3) {
      r = 0;
      g = 1. * (4 * igl::PI / 3 - angle) / (igl::PI / 3);
      b = 1;
    } else if (igl::PI * 4 / 3 <= angle && angle < igl::PI * 5 / 3) {
      r = 1. * (angle - 4 * igl::PI / 3) / (igl::PI / 3);
      g = 0;
      b = 1;
    } else if (igl::PI * 5 / 3 <= angle && angle <= 2 * igl::PI) {
      r = 1;
      g = 0;
      b = 1. * (2 * igl::PI - angle) / (igl::PI / 3);
    } else {
      r = 0, g = 0, b = 0;
    }

    assert(r <= 1 && r >= 0 && g <= 1 && g >= 0 && b <= 1 && b >= 0);
  }

  void angle_color(const Eigen::VectorXd &data, Eigen::MatrixX3d &color) {
    color.resize(data.rows(), 3);

    for (int i = 0; i < data.rows(); i++) {
      double angle = rp::round_pi(data(i));
      if (angle < 0) {
        angle += 2 * igl::PI;
      }

      assert(angle >= 0 && angle <= 2 * igl::PI);

      double r, g, b;
      get_rgb_from_angle(angle, r, g, b);
      color.row(i) << r, g, b;
    }
  }
}
