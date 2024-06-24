#include "./src/view/viewer.h"

int main(int argc, char *argv[])
{
  std::string path;

  if (argc < 2) {
    path = "";
  } else {
    path = argv[1];
  }
  rp::viewer::launch(path);
}
