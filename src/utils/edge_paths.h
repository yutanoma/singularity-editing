#pragma once

namespace rp {
  struct EdgePath {
    int edgeId;
    double ratio;

    EdgePath(const int &_edgeId, const double &_ratio);
  };
}
