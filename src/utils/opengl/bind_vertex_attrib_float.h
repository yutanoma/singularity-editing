#pragma once

#include <igl/opengl/gl.h>
#include <Eigen/Core>
#include <string>
namespace rp
{
  namespace opengl
  {
    // Bind a per-vertex array attribute and refresh its contents from an Eigen
    // matrix
    //
    // Inputs:
    //   program_shader  id of shader program
    //   name  name of attribute in vertex shader
    //   bufferID  id of buffer to bind to
    //   M  #V by dim matrix of per-vertex data
    //   refresh  whether to actually call glBufferData or just bind the buffer
    // Returns id of named attribute in shader
    GLint bind_vertex_attrib_float(
      const GLuint program_shader,
      const std::string &name, 
      GLuint bufferID, 
      const Eigen::VectorXf &M, 
      bool refresh);
  }
}
