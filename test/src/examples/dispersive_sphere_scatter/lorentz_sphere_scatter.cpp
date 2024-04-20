#include <memory>

#include "dispersive_sphere_scatter.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/material/dispersive_material.h"

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

  auto omega_p = xfdtd::Array1D<xfdtd::Real>{2 * xfdtd::constant::PI * 2e9};
  auto gamma = xfdtd::Array1D<xfdtd::Real>{xfdtd::constant::PI * 2e9};
  auto epsilon_inf = 2;
  auto epsilon_static = xfdtd::Array1D<xfdtd::Real>{5};
  auto nv = xfdtd::Array1D<xfdtd::Real>{xfdtd::constant::PI * 2e9};

  testCase(xfdtd::MLorentzMaterial::makeMLorentz(
               "lorentz_medium", epsilon_inf,
               (epsilon_static - epsilon_inf) * omega_p * omega_p, {0},
               omega_p * omega_p, 2 * nv, {1}),
           id);
}
