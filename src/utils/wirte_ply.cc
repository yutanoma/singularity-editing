#include "./write_ply.h"

#include <fstream>
#include <iostream>

namespace rp {
  void write_ply(const std::string &filename, const Eigen::MatrixX3i &F,
                 const Eigen::MatrixX3d &V, const Eigen::MatrixX3d &C) {
    std::ofstream writingFile;

    const int bufsize = 100000;
    char buf[bufsize];
    writingFile.rdbuf()->pubsetbuf(buf, bufsize);
    writingFile.open(filename, std::ios::trunc);

    if (!writingFile.is_open()) {
      std::cout << "file open failed" << std::endl;
      writingFile.close();
      return;
    }

    std::string text;

    writingFile << "ply" << std::endl;
    writingFile << "format ascii 1.0" << std::endl;
    writingFile << "element vertex " << V.rows() << std::endl;
    writingFile << "property float x" << std::endl;
    writingFile << "property float y" << std::endl;
    writingFile << "property float z" << std::endl;
    writingFile << "property uchar red" << std::endl;
    writingFile << "property uchar green" << std::endl;
    writingFile << "property uchar blue" << std::endl;
    writingFile << "element face " << F.rows() << std::endl;
    writingFile << "property list uchar int vertex_index" << std::endl;
    writingFile << "end_header" << std::endl;

    for (int i = 0; i < V.rows(); i++) {
      for (int j = 0; j < 3; j++) {
        writingFile << V(i, j) << " ";
      }
      for (int j = 0; j < 3; j++) {
        int v = (int)std::round(C(i, j) * 255);
        if (v < 0) {
          v = 0;
        } else if (v > 255) {
          v = 255;
        }
        writingFile << v;

        if (j != 2) {
          writingFile << " ";
        }
      }
      writingFile << std::endl;
    }

    for (int i = 0; i < F.rows(); i++) {
      writingFile << "3 ";
      for (int j = 0; j < 3; j++) {
        writingFile << F(i, j);

        if (j != 2) {
          writingFile << " ";
        }
      }
      writingFile << std::endl;
    }
  }
}
