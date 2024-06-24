#include <igl/remove_duplicate_vertices.h>
#include <igl/readSTL.h>
#include <igl/readOBJ.h>
#include <igl/unproject_onto_mesh.h>
#include <imgui/imgui.h>
#include <igl/PI.h>

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/meshio.h"

#include "../utils/opengl/glfw/Viewer.h"
#include "../utils/opengl/glfw/imgui/ImGuiMenu.h"
#include "../utils/opengl/glfw/imgui/ImGuiHelpers.h"

#include "../store/geometry.h"

#include "../action/types.h"
#include "../action/index_exchange.h"
#include "../action/edge_path.h"
#include "../action/parameterize.h"
#include "../action/vector_field.h"
#include "../action/rescaling.h"

#include "../../data/path.h"

#include "./set_geometry.h"
#include "./set_color.h"
#include "./set_angle.h"
#include "./set_path.h"
#include "./window.h"
#include "./mouse_down.h"
#include "./set_vector_field.h"
#include "./set_indices.h"

#include <future>
#include <fstream>
#include <iostream>
#include <chrono>

namespace rp {
namespace viewer {
bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

void launch(std::string filename = "")
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  // Plot the mesh
  rp::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  rp::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  Eigen::MatrixX3d V;
  Eigen::MatrixX3i F;
  Eigen::MatrixXi EV, EF, FE;

  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;

  // Load a mesh in OFF format
  // igl::readSTL(DATA_PATH "/models/bunny_mq_opened.stl", V, F, Eigen::MatrixX3d());
  std::string path = filename;
  if (path == "") {
    path = DATA_PATH "/models/bunny_mq.stl";
  }

  auto N = Eigen::MatrixX3d();
  if (ends_with(path, ".stl")) {
    auto res = igl::readSTL(path, V, F, N);
    std::tie(mesh, geometry) = geometrycentral::surface::readManifoldSurfaceMesh(path);

    assert(res);
  } else if (ends_with(path, ".obj")) {
    auto res = igl::readOBJ(path, V, F);
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(path);

    assert(res);
  } else {
    std::cerr << "Unsupported file format" << std::endl;
    return;
  }
  std::cout << "#faces=" << F.rows() << ", " << mesh->nFaces() << std::endl;

  geometry->requireVertexTangentBasis();
  auto vField = geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  Eigen::MatrixXd vf(vField.size(), 3);
  for (int i = 0; i < vField.size(); i++) {
    Vector3 basisX = geometry->vertexTangentBasis[i][0];
    Vector3 basisY = geometry->vertexTangentBasis[i][1];

    Vector2 field = vField[i];

    Vector3 v = field[0] * basisX + field[1] * basisY;

    vf.row(i) << v.x, v.y, v.z;
  }

  auto el = geometry->edgeLengths;
  double meanel = 0;
  for (int i = 0; i < el.size(); i++) {
    meanel += el[i];
  }
  meanel /= el.size();
  std::cout << "mean edge length: " << meanel << std::endl;

  std::cout << "initialization start" << std::endl;

  rp::store::initialize(F, V);
  rp::store::set_velocity(1 / meanel);
  rp::store::set_vectorfield(vf);
  rp::store::set_default_stripe_patterns();

  std::cout << "initialization end" << std::endl;

  auto pointSetter = rp::viewer::use_point_setter(viewer);
  int viewerId = 0;
  auto geometrySetter = rp::viewer::use_geometry_setter(viewer, V, F, viewerId);
  auto indicesSetter = rp::viewer::use_indices_setter(viewer, meanel / 3);
  auto pathSetter = rp::viewer::use_path_setter(viewer);
  auto fieldSetter = rp::viewer::use_vector_field_setter(viewer);
  auto fieldClearer = rp::viewer::use_vector_field_clearer(viewer);
  auto angleSetter = rp::viewer::use_angle_setter(viewer);

  auto _F = rp::store::get_f();
  auto _V = rp::store::get_v();
  geometrySetter(_F, _V);

  auto renderVectorfield = rp::action::use_show_vector_field(fieldSetter);
  auto parameterize = rp::action::use_parameterize(angleSetter, indicesSetter);
  auto indexExchange = rp::action::use_index_exchange_setter(angleSetter, indicesSetter);
  auto indexAnnihilate = rp::action::use_index_annihilation_setter(angleSetter, indicesSetter);
  auto addSingularityPair = rp::action::use_add_singularity_pair_setter(angleSetter, indicesSetter);
  auto computeRescaling = rp::action::use_rescaling_setter(angleSetter);
  auto clearRescaling = rp::action::use_rescaling_resetter(angleSetter);
  auto setBrushpath = rp::action::use_brush_path_setter(parameterize, pathSetter);
  auto addBrushPath = rp::action::use_brush_stroke_setter(pathSetter);

  parameterize();

  rp::viewer::SelectMode selectMode = rp::viewer::SelectMode::exchangeIndices;
  bool showVectorfield = false;

  auto resetInteraction = rp::viewer::use_reset_interaction();

  auto customWindowCallback = rp::viewer::use_custom_window(viewer, menu, selectMode, showVectorfield, resetInteraction, renderVectorfield, fieldClearer, parameterize, computeRescaling, clearRescaling);
  menu.callback_draw_custom_window = customWindowCallback;

  auto mouseDownCallback = use_mouse_down(V, F, selectMode, geometrySetter, pathSetter, pointSetter, setBrushpath, setBrushpath, indexExchange, indexAnnihilate, addSingularityPair, addBrushPath);
  viewer.callback_mouse_down = mouseDownCallback;

  viewer.data().set_face_based(true);
  viewer.core().background_color << 255, 255, 255, 1;

  viewer.launch_init(true, false);
  viewer.launch_rendering(true);
  viewer.launch_shut();
}
}
}
