#include "flatten_vect.h"

namespace rp {
template <typename T>
T flatten_vect(std::vector<T> vect) {
  int new_size = 0;
  int cols = 3;

  for (int i = 0; i < vect.size(); i++) {
    new_size += vect[i].rows();
    cols = vect[i].cols();
  }

  T result;
  result.resize(new_size, cols);

  int count = 0;
  for (int i = 0; i < vect.size(); i++) {
    for (int j = 0; j < vect[i].rows(); j++) {
      for (int k = 0; k < cols; k++) {
        result(count, k) = vect[i](j, k);
      }

      count++;
    }
  }

  return result;
};

// 明示的なインスタンス
template Eigen::MatrixX3d flatten_vect<Eigen::MatrixX3d>(
    std::vector<Eigen::MatrixX3d> vect);

template Eigen::MatrixX3i flatten_vect<Eigen::MatrixX3i>(
    std::vector<Eigen::MatrixX3i> vect);
}  // namespace rp
