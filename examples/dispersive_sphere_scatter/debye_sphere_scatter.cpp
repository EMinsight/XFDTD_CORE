#include "dispersive_sphere_scatter.h"
#include "xfdtd/util/constant.h"

int main(int argc, char *argv[]) {
  int id = 1;
  if (argc != 2) {
    id = 1;
  }

  if (argc == 2) {
    id = std::stoi(argv[1]);
    if (id < 0 || 2 <= id) {
      id = 1;
    }
  }

  testCase(std::make_unique<xfdtd::DebyeMedium>(xfdtd::DebyeMedium{
               "debye_medium", 2, {7}, {2e-9 / (2 * xfdtd::constant::PI)}}),
           id);
}
