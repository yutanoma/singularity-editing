#include "./triangle_angles.h"

namespace rp {
  void get_angle(const double &a, const double &b, const double &c,
                 double &A, double &B, double &C) {
    double s = (a + b + c) / 2;
    double area = std::sqrt(s * (s - a) * (s - b) * (s - c));

    double h = 2 * area / a;
    B = std::asin(h / c);
    C = std::asin(h / b);
    A = igl::PI - B - C;

    assert(std::min({A, B, C}) >= 0);
  }

  void triangle_angles(const double &a, const double &b, const double &c,
                       double &A, double &B, double &C) {
    double max = std::max({a, b, c});
    double min = std::min({a, b, c});
    assert(min >= 0);

    if (a == max) {
      get_angle(a, b, c, A, B, C);
    } else if (b == max) {
      get_angle(b, a, c, B, A, C);
    } else if (c == max) {
      get_angle(c, a, b, C, A, B);
    } else {
      assert(false && "no max in abc");
    }
  }
}
