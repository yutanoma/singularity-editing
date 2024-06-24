// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "../Viewer.h"
#include "../ViewerPlugin.h"
#include <igl/igl_inline.h>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

// Forward declarations
struct ImGuiContext;

namespace rp
{
namespace opengl
{
namespace glfw
{
namespace imgui
{

class ImGuiMenu : public rp::opengl::glfw::ViewerPlugin
{
protected:
  // Hidpi scaling to be used for text rendering.
  float hidpi_scaling_;

  // Ratio between the framebuffer size and the window size.
  // May be different from the hipdi scaling!
  float pixel_ratio_;

  // ImGui Context
  ImGuiContext * context_ = nullptr;

public:
  virtual void init(rp::opengl::glfw::Viewer *_viewer) override;

  virtual void reload_font(int font_size = 13);

  virtual void shutdown() override;

  virtual bool pre_draw() override;

   virtual bool post_draw() override;

  virtual void post_resize(int width, int height) override;

  // Mouse IO
  virtual bool mouse_down(int button, int modifier) override;

  virtual bool mouse_up(int button, int modifier) override;

  virtual bool mouse_move(int mouse_x, int mouse_y) override;

  virtual bool mouse_scroll(float delta_y) override;

  // Keyboard IO
  virtual bool key_pressed(unsigned int key, int modifiers) override;

  virtual bool key_down(int key, int modifiers) override;

  virtual bool key_up(int key, int modifiers) override;

  // Draw menu
  virtual void draw_menu();

  // Can be overwritten by `callback_draw_viewer_window`
  virtual void draw_viewer_window();

  // Can be overwritten by `callback_draw_viewer_menu`
  virtual void draw_viewer_menu();

  // Can be overwritten by `callback_draw_custom_window`
  virtual void draw_custom_window() { }

  // Easy-to-customize callbacks
  std::function<void(void)> callback_draw_viewer_window;
  std::function<void(void)> callback_draw_viewer_menu;
  std::function<void(void)> callback_draw_custom_window;

  void draw_text(
    Eigen::Vector3d pos,
    Eigen::Vector3d normal,
    const std::string &text,
    const Eigen::Vector4f color = Eigen::Vector4f(0,0,0.04,1)); // old default color

  float pixel_ratio();

  float hidpi_scaling();

  float menu_scaling() { return hidpi_scaling_ / pixel_ratio_; }
};

} // end namespace
} // end namespace
} // end namespace
} // end namespace
