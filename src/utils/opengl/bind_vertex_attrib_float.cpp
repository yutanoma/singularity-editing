#include "bind_vertex_attrib_float.h"

GLint rp::opengl::bind_vertex_attrib_float(
  const GLuint program_shader,
  const std::string &name, 
  GLuint bufferID,
  const Eigen::VectorXf &M, 
  bool refresh)
{
  GLint id = glGetAttribLocation(program_shader, name.c_str());
  if (id < 0)
    return id;
  if (M.size() == 0)
  {
    glDisableVertexAttribArray(id);
    return id;
  }
  glBindBuffer(GL_ARRAY_BUFFER, bufferID);
  if (refresh)
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*M.size(), M.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(id, M.cols(), GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(id);
  return id;
}
