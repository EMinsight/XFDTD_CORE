
#include <xfdtd/common/constant.h>

#include "dispersive_sphere_scatter.h"

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

  testCase(std::make_unique<xfdtd::DrudeMedium>(xfdtd::DrudeMedium{
               "drude_medium", 4, {1.5e9 * 2 * xfdtd::constant::PI}, {1e9}}),
           id);
}
