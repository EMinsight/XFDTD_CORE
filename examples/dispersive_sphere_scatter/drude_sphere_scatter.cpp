#include <xfdtd/common/constant.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/material/dispersive_material.h>

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

  xfdtd::Array1D<xfdtd::Real> omega_p = {xfdtd::constant::PI * 1e9};
  xfdtd::Array1D<xfdtd::Real> gamma = {xfdtd::constant::PI * 1.2e9};

  testCase(xfdtd::LinearDispersiveMaterial::makeMLorentz(
               "drude_m_lor", 4, {omega_p * omega_p}, {0}, {0}, {gamma}, {1}),
           id);

  // testCase(
  //     xfdtd::DrudeMedium::makeDrudeMedium("drude_medium", 4, omega_p, gamma),
  //     id);
}
