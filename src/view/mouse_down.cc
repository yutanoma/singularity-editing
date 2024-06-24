#include "./mouse_down.h"

#include <igl/unproject_onto_mesh.h>

namespace rp {
namespace viewer {
namespace {
  int prevVid = -1;
  int prevVid2 = -1;
  int prevPrevVid = -1;
  int prevFid = -1;
  bool mouseDown = false;
  SelectMode lastSelected = SelectMode::undefined;
};

bool exchange_indices(const int &fid, rp::action::IndexExchangeStateSetter &exchangeIndexSetter) {
  if (prevFid == -1 || lastSelected != SelectMode::exchangeIndices) {
    prevFid = fid;
  } else {
    exchangeIndexSetter(fid, prevFid);
    prevFid = -1;
  }

  return true;
}

bool annihilate_indices(const int &fid, rp::action::IndexExchangeStateSetter &exchangeIndexSetter) {
  if (prevFid == -1 || lastSelected != SelectMode::annihilateIndices) {
    prevFid = fid;
  } else {
    exchangeIndexSetter(fid, prevFid);
    prevFid = -1;
  }

  return true;
}

bool add_singularities(const int &fid, rp::action::IndexExchangeStateSetter &singularityPairSetter) {
  std::cout << fid << ", " << prevFid << std::endl;
  if (prevFid == -1 || lastSelected != SelectMode::addSingularities) {
    prevFid = fid;
  } else {
    singularityPairSetter(fid, prevFid);
    prevFid = -1;
  }

  return true;
}

bool draw_geodesic_brush_path(const int &vid, rp::action::EdgePathStateSetter &setBrushpath) {
  if (prevVid == -1 || lastSelected != SelectMode::brushPath) {
    prevVid = vid;
  } else {
    setBrushpath(prevVid, vid);
    prevVid = vid;
  }

  return true;
}

bool draw_brush_path(const int &vid, const int &secondVid, const double &ratio, rp::action::EdgePathBrushStateSetter &setBrushpath) {
  int v0 = std::min(vid, secondVid), v1 = std::max(vid, secondVid);

  bool isSomeEqual = v0 == prevVid || v0 == prevVid2 || v1 == prevVid || v1 == prevVid2;
  bool isNotSame = std::min(prevVid, prevVid2) != v0 && std::max(prevVid, prevVid2) != v1;

  if ((isSomeEqual && isNotSame) || (prevVid == -1 && prevVid2 == -1)) {
    setBrushpath(v0, v1, ratio);
    prevVid = v0;
    prevVid2 = v1;
  }

  return true;
}

bool draw_prescription_path() {
  return false;
}

bool default_operation() {
  prevVid = -1;
  prevVid2 = -1;
  prevPrevVid = -1;
  prevFid = -1;
  lastSelected = SelectMode::undefined;
  mouseDown = false;
  return false;
}

MouseDownCallback use_mouse_down(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                                 SelectMode &selectMode,
                                 rp::action::GeometrySetter &geometrySetter,
                                 rp::action::PathSetter &pathSetter,
                                 rp::action::PointsSetter &pointSetter,
                                 rp::action::EdgePathStateSetter &setGeodesicBrushPath,
                                 rp::action::EdgePathStateSetter &setGeodesicPrescriptionPath,
                                 rp::action::IndexExchangeStateSetter &exchangeIndices,
                                 rp::action::IndexExchangeStateSetter &annihilateIndices,
                                 rp::action::IndexExchangeStateSetter &addSingularityPair,
                                 rp::action::EdgePathBrushStateSetter &setBrushpath) {
  auto mouseDownCallback = [&](rp::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;

    if (igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport,
                                 V, F, fid, bc)) {
      int max, min, secondMax;
      bc.maxCoeff(&max);
      bc.minCoeff(&min);
      secondMax = 3 - max - min;
      if (secondMax > 2 || secondMax < 0) {
        assert(min == max);
        secondMax = (min + 1) % 3;
      }
      double ratio = bc(secondMax) / (bc(max) + bc(secondMax));
      int vid = F(fid,max);
      int secondVid = F(fid, secondMax);
      bool result = false;

      std::cout << "fid: " << fid << ", vid: " << vid << ", secondvid: " << secondVid << std::endl;

      switch (selectMode) {
        case SelectMode::exchangeIndices:
          result = exchange_indices(fid, exchangeIndices);
          break;
        case SelectMode::annihilateIndices:
          result = annihilate_indices(fid, annihilateIndices);
          break;
        case SelectMode::addSingularities:
          result = add_singularities(fid, addSingularityPair);
          break;
        case SelectMode::brushPath:
          result = draw_geodesic_brush_path(vid, setGeodesicBrushPath);
          break;
        case SelectMode::indexPrescriptionPath:
          result = draw_prescription_path();
          break;
        default:
          return default_operation();
      }

      lastSelected = selectMode;
      return result;
    } else {
      return default_operation();
    }
  };

  return mouseDownCallback;
}

std::function<void()> use_reset_interaction() {
  std::function<void()> resetInteraction = [&]()
  {
    default_operation();
  };

  return resetInteraction;
}
}
}