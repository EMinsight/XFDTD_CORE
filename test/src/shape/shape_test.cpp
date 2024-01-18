#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/cylinder.h"
#include "xfdtd/shape/sphere.h"
void cubeTest() {
  auto cube_0{xfdtd::Cube{xfdtd::Vector{0, 0, 0}, xfdtd::Vector{1, 1, 1}}};
  assert(cube_0.isInside(0.5, 0.5, 0.5));
  assert(!cube_0.isInside(1.5, 0.5, 0.5));
  assert(!cube_0.isInside(0.5, 1.5, 0.5));
  assert(!cube_0.isInside(0.5, 0.5, 1.5));
  assert(!cube_0.isInside(-0.5, 0.5, 0.5));
  assert(!cube_0.isInside(0.5, -0.5, 0.5));
}

void cylinderTest() {
  auto cylinder_0{
      xfdtd::Cylinder{xfdtd::Vector{0.5, 0, 0}, 1, 1, xfdtd::Axis::XP}};
  assert(cylinder_0.isInside(0.5, 0.5, 0.5));
  assert(!cylinder_0.isInside(1.5, 0.5, 0.5));
  assert(!cylinder_0.isInside(0.5, 1.5, 0.5));
  assert(!cylinder_0.isInside(0.5, 0.5, 1.5));
  assert(!cylinder_0.isInside(-0.5, 0.5, 0.5));
  assert(cylinder_0.isInside(0.5, -0.5, 0.5));
}

void sphereTest() {
  auto sphere_0{xfdtd::Sphere{xfdtd::Vector{0, 0, 0}, 1}};
  assert(sphere_0.isInside(0.5, 0.5, 0.5));
  assert(!sphere_0.isInside(1.5, 0.5, 0.5));
  assert(!sphere_0.isInside(0.5, 1.5, 0.5));
  assert(!sphere_0.isInside(0.5, 0.5, 1.5));
  assert(sphere_0.isInside(-0.5, 0.5, 0.5));
  assert(sphere_0.isInside(0.5, -0.5, 0.5));
}

int main() {
  cubeTest();
  cylinderTest();
  sphereTest();
  return 0;
}
