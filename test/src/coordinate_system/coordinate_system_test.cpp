#include <xfdtd/common/constant.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

#include <cassert>

#include "util/float_compare.h"

void vectorTest() {
  xfdtd::Vector vector{1, 2, 3};
  assert(xfdtd::floatCompareEqual(1, vector.x()));
  assert(xfdtd::floatCompareEqual(2, vector.y()));
  assert(xfdtd::floatCompareEqual(3, vector.z()));
  assert(xfdtd::floatCompareEqual(std::sqrt(14), vector.normL2()));

  xfdtd::Vector normalized_vector{xfdtd::Vector::normalized(vector)};
  assert(xfdtd::floatCompareEqual(1 / std::sqrt(14), normalized_vector.x()));
  assert(xfdtd::floatCompareEqual(2 / std::sqrt(14), normalized_vector.y()));
  assert(xfdtd::floatCompareEqual(3 / std::sqrt(14), normalized_vector.z()));
  assert(xfdtd::floatCompareEqual(1, normalized_vector.normL2()));

  xfdtd::Vector spherical_vector{xfdtd::Vector::fromSpherical(
      2, xfdtd::constant::PI / 4, xfdtd::constant::PI / 2)};
  assert(xfdtd::floatCompareEqual(
      2 * std::sin(xfdtd::constant::PI / 4) * std::cos(xfdtd::constant::PI / 2),
      spherical_vector.x()));
  assert(xfdtd::floatCompareEqual(
      2 * std::sin(xfdtd::constant::PI / 4) * std::sin(xfdtd::constant::PI / 2),
      spherical_vector.y()));
  assert(xfdtd::floatCompareEqual(2 * std::cos(xfdtd::constant::PI / 4),
                                  spherical_vector.z()));

  auto cylindrical_vector{
      xfdtd::Vector::fromCylindrical(2, xfdtd::constant::PI / 4, 2)};
  assert(xfdtd::floatCompareEqual(2 * std::cos(xfdtd::constant::PI / 4),
                                  cylindrical_vector.x()));
  assert(xfdtd::floatCompareEqual(2 * std::sin(xfdtd::constant::PI / 4),
                                  cylindrical_vector.y()));
  assert(xfdtd::floatCompareEqual(2, cylindrical_vector.z()));
  assert(xfdtd::floatCompareEqual(std::sqrt(8), cylindrical_vector.normL2()));

  // compare test
  auto vector1{xfdtd::Vector{xfdtd::constant::PI / 2, xfdtd::constant::PI,
                             xfdtd::constant::PI * 2}};
  auto vector2{xfdtd::Vector{xfdtd::constant::PI / 2, xfdtd::constant::PI,
                             xfdtd::constant::PI * 2}};
  assert(vector1 == vector2);
  assert(!(vector1 != vector2));
  auto vector3{xfdtd::Vector{xfdtd::constant::PI / 2, xfdtd::constant::PI,
                             xfdtd::constant::PI * 2 - 0.1}};
  assert(vector1 != vector3);
  assert(!(vector1 == vector3));
}

void axisTest() {
  auto vector{xfdtd::Vector{1, 2, 3}};
  auto axis{xfdtd::Axis{vector}};

  assert(xfdtd::floatCompareEqual(axis.normL2(), 1));
  assert(!(axis == xfdtd::Axis::XN));

  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::XP) ==
         xfdtd::Axis::XP);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::XN) ==
         xfdtd::Axis::XN);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::YP) ==
         xfdtd::Axis::YP);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::YN) ==
         xfdtd::Axis::YN);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::ZP) ==
         xfdtd::Axis::ZP);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::ZN) ==
         xfdtd::Axis::ZN);

  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::XP) !=
         xfdtd::Axis::XN);
  assert(xfdtd::Axis::fromDirectionToAxis(xfdtd::Axis::Direction::XN) !=
         xfdtd::Axis::XP);
  auto axis_1{xfdtd::Axis{vector}};
  assert(axis_1 == axis);

  auto xp{xfdtd::Axis{xfdtd::Vector{1, 0, 0}}};
  assert(xp == xfdtd::Axis::XP);
  assert(xp != xfdtd::Axis::XN);

  assert(xp == xfdtd::Axis::XYZ::X);
  assert(xp != xfdtd::Axis::XYZ::Y);
  assert(xp != xfdtd::Axis::XYZ::Z);

  bool exception_thrown{false};
  try {
    auto axis_err{xfdtd::Axis{xfdtd::Vector{0, 0, 0}}};
    exception_thrown = false;
  } catch (std::exception& e) {
    exception_thrown = true;
  }
  assert(exception_thrown);
}

int main() {
  vectorTest();
  axisTest();
  return 0;
}
