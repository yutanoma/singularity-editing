#include "./edge_path.h"

#include "../dispatcher/edge_path.h"
#include <iostream>

namespace rp {
namespace action {
EdgePathStateSetter use_brush_path_setter(ParameterizationStateSetter &parameterize, PathSetter &pathSetter) {
  auto edgePathSetter = [&](const int &prevVid, const int &nextVid)
  {
    std::cout << prevVid << ", " << nextVid << std::endl;
    Eigen::MatrixX3d newPath;
    rp::dispatcher::set_brush_paths(prevVid, nextVid, newPath);
    Eigen::MatrixX3d C(1, 3);
    C << 1, .5, 0;
    pathSetter(newPath, C);
    // computeVectorfield();
    // parameterize();
  };

  return edgePathSetter;
}

EdgePathBrushStateSetter use_brush_stroke_setter(PathSetter &pathSetter) {
  auto strokeSetter = [&](const int &vid0, const int &vid1, const double &ratio) {
    Eigen::MatrixX3d newPath;
    rp::dispatcher::set_brush_stroke(vid0, vid1, ratio, newPath);
    Eigen::MatrixX3d C(1, 3);
    C << 1, .5, 0;
    pathSetter(newPath, C);
  };

  return strokeSetter;
}
}
}

