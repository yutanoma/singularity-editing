// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

// Coverts mesh data inside a igl::ViewerData class in an OpenGL
// compatible format The class includes a shader and the opengl calls to plot
// the data

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace rp
{
namespace opengl
{

class MeshGL
{
public:
  typedef unsigned int GLuint;
  typedef unsigned int GLint;

  enum DirtyFlags
  {
    DIRTY_NONE           = 0x0000,
    DIRTY_POSITION       = 0x0001,
    DIRTY_UV             = 0x0002,
    DIRTY_NORMAL         = 0x0004,
    DIRTY_AMBIENT        = 0x0008,
    DIRTY_DIFFUSE        = 0x0010,
    DIRTY_SPECULAR       = 0x0020,
    DIRTY_TEXTURE        = 0x0040,
    DIRTY_FACE           = 0x0080,
    DIRTY_MESH           = 0x00FF,
    DIRTY_OVERLAY_LINES  = 0x0100,
    DIRTY_OVERLAY_POINTS = 0x0200,
    DIRTY_VERTEX_LABELS  = 0x0400,
    DIRTY_FACE_LABELS    = 0x0800,
    DIRTY_CUSTOM_LABELS  = 0x1000,
    DIRTY_ALL            = 0xFFFF
  };

  bool is_initialized = false;
  GLuint vao_mesh;
  GLuint vao_overlay_lines;
  GLuint vao_overlay_points;
  GLuint shader_mesh;
  GLuint shader_overlay_lines;
  GLuint shader_overlay_points;
  GLuint shader_text;

  GLuint vbo_V; // Vertices of the current mesh (#V x 3)
  GLuint vbo_V_uv; // UV coordinates for the current mesh (#V x 2)
  GLuint vbo_V_normals; // Vertices of the current mesh (#V x 3)
  GLuint vbo_V_ambient; // Ambient material  (#V x 3)
  GLuint vbo_V_diffuse; // Diffuse material  (#V x 3)
  GLuint vbo_V_specular; // Specular material  (#V x 3)
  GLuint vbo_V_angle; // 色を割り当てるときの角度値
  GLuint vbo_V_index; // 特異点の指数
  GLuint vbo_V_indexInFace; // その点の属する面の中で何番目の点か

  GLuint vbo_F; // Faces of the mesh (#F x 3)
  GLuint vbo_tex; // Texture

  GLuint vbo_lines_F;         // Indices of the line overlay
  GLuint vbo_lines_V;         // Vertices of the line overlay
  GLuint vbo_lines_V_colors;  // Color values of the line overlay
  GLuint vbo_points_F;        // Indices of the point overlay
  GLuint vbo_points_V;        // Vertices of the point overlay
  GLuint vbo_points_V_colors; // Color values of the point overlay

  // Temporary copy of the content of each VBO
  typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXf;
  RowMatrixXf V_vbo;
  RowMatrixXf V_normals_vbo;
  RowMatrixXf V_ambient_vbo;
  RowMatrixXf V_diffuse_vbo;
  RowMatrixXf V_specular_vbo;
  RowMatrixXf V_uv_vbo;
  Eigen::VectorXf V_angle_vbo;
  Eigen::VectorXf V_index_vbo;
  Eigen::VectorXf V_indexInFace_vbo;
  bool enable_angle_color;
  RowMatrixXf lines_V_vbo;
  RowMatrixXf lines_V_colors_vbo;
  RowMatrixXf points_V_vbo;
  RowMatrixXf points_V_colors_vbo;

  // Text Rendering
  struct TextGL
  { 
    uint32_t dirty_flag;
    GLuint vao_labels;
    GLuint vbo_labels_pos;
    GLuint vbo_labels_characters;
    GLuint vbo_labels_offset; 
    GLuint vbo_labels_indices;
    RowMatrixXf label_pos_vbo;
    RowMatrixXf label_char_vbo;
    RowMatrixXf label_offset_vbo;
    Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> label_indices_vbo;
    void init_buffers();
    void free_buffers();
  };
  TextGL vertex_labels;
  TextGL face_labels;  
  TextGL custom_labels;
  GLuint font_atlas;

  int tex_u;
  int tex_v;
  GLint tex_filter;
  GLint tex_wrap;
  Eigen::Matrix<char,Eigen::Dynamic,1> tex;

  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> F_vbo;
  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> lines_F_vbo;
  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> points_F_vbo;

  // Marks dirty buffers that need to be uploaded to OpenGL
  uint32_t dirty;

  MeshGL();

  // Initialize shaders and buffers
  void init();

  // Release all resources
  void free();

  // Create a new set of OpenGL buffer objects
  void init_buffers();

  // Bind the underlying OpenGL buffer objects for subsequent mesh draw calls
  void bind_mesh();

  /// Draw the currently buffered mesh (either solid or wireframe)
  void draw_mesh(bool solid);

  // Bind the underlying OpenGL buffer objects for subsequent line overlay draw calls
  void bind_overlay_lines();

  /// Draw the currently buffered line overlay
  void draw_overlay_lines();

  // Bind the underlying OpenGL buffer objects for subsequent point overlay draw calls
  void bind_overlay_points();

  /// Draw the currently buffered point overlay
  void draw_overlay_points();

  // Text Binding and Draw functions
  void init_text_rendering();
  void bind_labels(const TextGL& labels);
  void draw_labels(const TextGL& labels);

  // Release the OpenGL buffer objects
  void free_buffers();

};

}
}
