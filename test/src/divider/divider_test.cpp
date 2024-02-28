#include "divider/divider.h"

#include <iostream>

int main() {
  auto problem =
      xfdtd::Divider::Task<std::size_t>{{0, 100}, {6, 51}, {0, 74}};
  auto num_procs = 4;
  auto res = xfdtd::Divider::divide(problem, num_procs);
  for (auto&& r : res) {
    std::cout << r.toString() << "\n";
  }
  res = xfdtd::Divider::divide(problem, num_procs, xfdtd::Divider::Type::X);
  for (auto&& r : res) {
    std::cout << r.toString() << "\n";
  }
  return 0;
}