#include <memory>

#include "dispersive_sphere_scatter.h"

int main(int argc, char* argv[]) {
  int id = 0;

  if (argc != 2) {
    id = 1;
  }

  if (argc == 2) {
    id = std::stoi(argv[1]);
    if (id < 0 || 2 <= id) {
      id = 1;
    }
  }

  testCase(std::make_shared<xfdtd::LorentzMedium>(
               xfdtd::LorentzMedium{"lorentz_medium",
                                    2,
                                    {5},
                                    {2 * xfdtd::constant::PI * 2e9},
                                    {xfdtd::constant::PI * 2e9}}),
           id);
}
