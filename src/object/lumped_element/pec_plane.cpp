#include <xfdtd/object/lumped_element/pec_plane.h>

#include <xtensor.hpp>

#include "xfdtd/object/object.h"

namespace xfdtd {

PecPlane::PecPlane(std::string name, std::unique_ptr<Cube> cube)
    : Object{std::move(name), std::move(cube), Material::createPec()} {}

void PecPlane::correctMaterialSpace() {
  auto grid_box{gridBoxPtr()};
  auto em_property{materialPtr()->emProperty()};
  auto eps{em_property.epsilon()};
  auto mu{em_property.mu()};
  auto sigma_e{em_property.sigmaE()};
  auto sigma_m{em_property.sigmaM()};
  auto calc_param{calculationParamPtr()};

  if (grid_box->origin().i() == grid_box->end().i()) {
    auto is{grid_box->origin().i()};
    auto js{grid_box->origin().j()};
    auto ks{grid_box->origin().k()};
    auto ie{grid_box->end().i()};
    auto je{grid_box->end().j()};
    auto ke{grid_box->end().k()};
    xt::view(calc_param->materialParam()->sigmaEY(), is, xt::range(js, je),
             xt::range(ks, ke + 1)) = sigma_e;

    xt::view(calc_param->materialParam()->sigmaEZ(), is, xt::range(js, je + 1),
             xt::range(ks, ke)) = sigma_e;
    return;
  }

  if (grid_box->origin().j() == grid_box->end().j()) {
    auto is{grid_box->origin().i()};
    auto js{grid_box->origin().j()};
    auto ks{grid_box->origin().k()};
    auto ie{grid_box->end().i()};
    auto je{grid_box->end().j()};
    auto ke{grid_box->end().k()};
    xt::view(calc_param->materialParam()->sigmaEZ(), xt::range(is, ie + 1), js,
             xt::range(ks, ke)) = sigma_e;

    xt::view(calc_param->materialParam()->sigmaEX(), xt::range(is, ie), js,
             xt::range(ks, ke + 1)) = sigma_e;
    return;
  }

  if (grid_box->origin().k() == grid_box->end().k()) {
    auto is{grid_box->origin().i()};
    auto js{grid_box->origin().j()};
    auto ks{grid_box->origin().k()};
    auto ie{grid_box->end().i()};
    auto je{grid_box->end().j()};
    auto ke{grid_box->end().k()};
    xt::view(calc_param->materialParam()->sigmaEX(), xt::range(is, ie),
             xt::range(js, je + 1), ks) = sigma_e;

    xt::view(calc_param->materialParam()->sigmaEY(), xt::range(is, ie + 1),
             xt::range(js, je), ks) = sigma_e;
    return;
  }

  throw XFDTDObjectException(name() + " PecPlane is not a plane!");
}

}  // namespace xfdtd
