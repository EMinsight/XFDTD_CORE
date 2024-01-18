#include <iostream>
#include <xtensor.hpp>

int main() {
  xt::xarray<double> m = xt::random::randn<double>({4, 4, 4});
  auto v{xt::view(m, xt::range(1, 4), xt::all(), xt::all())};
  std::cout << "m: " << m << "\n";
  std::cout << "v: " << v << "\n";

  v *= 2.0;
  std::cout << "do v *= 2.0"
            << "\n";
  std::cout << "m: " << m << "\n";
  std::cout << "v: " << v << "\n";
  std::cout << "\n";

  // random
  auto k{xt::xarray<double>{2, 2,2}};
  std::cout << "k: " << k << "\n";

  // do add
  std::cout << v / xt::view(k, xt::all(), xt::newaxis(), xt::newaxis()) << "\n";
  v /= xt::view(k, xt::all(), xt::newaxis(), xt::newaxis());
  std::cout << "M: \n" <<  xt::view(m, xt::range(1, 4), xt::all(), xt::all()) << "\n";
}
