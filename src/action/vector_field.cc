#include "./vector_field.h"

#include <iostream>

namespace rp {
namespace action {
  VectorFieldStateSetter use_show_vector_field(FieldSetter &fieldSetter) {
    VectorFieldStateSetter showVectorField = [&]()
    {
      try
      {
        auto vf = rp::store::get_vector_field();

        auto F = rp::store::get_f();
        auto V = rp::store::get_v();

        Eigen::MatrixXd vfOnFace(F.rows(), 3);
        vfOnFace.setZero();

        for (int i = 0; i < F.rows(); i++) {
          for (int j = 0; j < 3; j++) {
            int vid = F(i, j);
            vfOnFace.row(i) = vfOnFace.row(i) + vf.row(vid);
          }
          vfOnFace.row(i) = vfOnFace.row(i) / 3;
        }

        fieldSetter(vfOnFace, F, V);
      }
      catch (std::exception &e)
      {
        std::cout << e.what() << std::endl;
      }
    };

    return showVectorField;
  }
}
}
