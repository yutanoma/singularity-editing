#pragma once
#include <igl/igl_inline.h>
#include "./ImGuiMenu.h"
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imguizmo/ImGuizmo.h>
#include <Eigen/Dense>
#include <vector>

namespace rp{ namespace opengl{ namespace glfw{ namespace imgui{

class SelectionPlugin: public rp::opengl::glfw::imgui::ImGuiMenu
{
public:
  // customizable hotkeys
  std::string MARQUEE_KEY = "Mm";
  // leave 'L' for show_lines in viewer
  std::string LASSO_KEY = "l";
  std::string OFF_KEY = "Vv";
  enum Mode
  {
    OFF                 = 0,
    RECTANGULAR_MARQUEE = 1,
    ELLIPTICAL_MARQUEE  = 2,
    POLYGONAL_LASSO     = 3,
    LASSO               = 4,
    NUM_MODES           = 5
  } mode = RECTANGULAR_MARQUEE;
  bool is_down = false;
  bool has_moved_since_down = false;
  bool is_drawing = false;
  // min and max corners of 2D rectangular marquee
  Eigen::Matrix<float,2,2> M = Eigen::Matrix<float,2,2>::Zero();
  // list of points of 2D lasso marquee
  std::vector<Eigen::RowVector2f> L;
  // callback called when slection is completed (usually on mouse_up)
  std::function<void(void)> callback;
  // callback called after mode is changed 
  std::function<void(Mode)> callback_post_mode_change;
  // whether rotating, translating or scaling
  ImGuizmo::OPERATION operation;
  // stored transformation
  Eigen::Matrix4f T;
  // Initilize with rotate operation on an identity transform (at origin)
  SelectionPlugin():operation(ImGuizmo::ROTATE),T(Eigen::Matrix4f::Identity()){};
  virtual void init(rp::opengl::glfw::Viewer *_viewer) override;
  virtual bool pre_draw() override;
  virtual bool post_draw() override;
  virtual bool mouse_down(int button, int modifier) override;
  virtual bool mouse_up(int button, int modifier) override;
  virtual bool mouse_move(int mouse_x, int mouse_y) override;
  virtual bool key_pressed(unsigned int key, int modifiers) override;
  void clear();
  // helpers
  static void circle(const Eigen::Matrix<float,2,2> & M,  std::vector<Eigen::RowVector2f> & L);
  static void rect(const Eigen::Matrix<float,2,2> & M,  std::vector<Eigen::RowVector2f> & L);
  static Eigen::RowVector2f xy(const Viewer * v);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}}}}
