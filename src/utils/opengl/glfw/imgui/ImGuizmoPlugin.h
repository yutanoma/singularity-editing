#pragma once
#include <igl/igl_inline.h>
#include "ImGuiMenu.h"
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imguizmo/ImGuizmo.h>
#include <Eigen/Dense>

namespace rp{ namespace opengl{ namespace glfw{ namespace imgui{

class ImGuizmoPlugin : public rp::opengl::glfw::imgui::ImGuiMenu
{
public:
  // callback(T) called when the stored transform T changes
  std::function<void(const Eigen::Matrix4f &)> callback;
  // Whether to display
  bool visible = true;
  // whether rotating, translating or scaling
  ImGuizmo::OPERATION operation;
  // stored transformation
  Eigen::Matrix4f T;
  // Initilize with rotate operation on an identity transform (at origin)
  ImGuizmoPlugin():operation(ImGuizmo::ROTATE),T(Eigen::Matrix4f::Identity()){};
  /////////////////////////////////////////////////////////////////////////////
  // Boilerplate
  virtual void init(rp::opengl::glfw::Viewer *_viewer) override;
  virtual bool pre_draw() override;
  /////////////////////////////////////////////////////////////////////////////
  virtual bool post_draw() override;
};

}}}}
