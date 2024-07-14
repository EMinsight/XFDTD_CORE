# XFDTD CORE

This is a part of the XFDTD project. It contains the core functionality of the XFDTD project.

## Feature

1. Support 1D, 2D and 3D simulation.
2. Support PML boundary.
3. Support TFSF source.
4. Support lumped element.
5. Support dispersive material.(Lorentz, Drude, Debye and Modified Lorentz)
6. Parallel computing support.(C++ Standard Thread, MPI and CUDA)
7. cross-platform support.(Linux, Windows and MacOS)

## Getting Stared

You will require the the following libraries on your system to compile:[xtl](https://github.com/xtensor-stack/xtl), [xtensor](https://github.com/xtensor-stack/xtensor).

C++20 is required to compile the library.

### Install from source

```bash
cmake -S . -B build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --target install
```

make `MPI` to `ON` if you want to use MPI.

```bash
cmake -DXFDTD_CORE_WITH_MPI=ON -S . -B build
```

### Use in your project

Here is a example to simulate a dielectric sphere scatter a plane wave.

Assuming you have installed the library in `/path/to/install`, you can use the following `CMakeLists.txt` to use the library in your project. The folder structure:

```tree
.
├── CMakeLists.txt
└── main.cpp
```

The `CMakeLists.txt` file:

```cmake
cmake_minimum_required(VERSION 3.20)

project(xfdtd_first_example VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(xfdtd_core REQUIRED)

set(XFDTD_EXAMPLE_LIBRARIES xfdtd::xfdtd_core)

# check for MPI
if(XFDTD_CORE_WITH_MPI)
  set(XFDTD_EXAMPLE_LIBRARIES ${XFDTD_EXAMPLE_LIBRARIES} MPI::MPI_CXX)
endif()

add_executable(xfdtd_first_example main.cpp)

target_link_libraries(xfdtd_first_example PRIVATE ${XFDTD_EXAMPLE_LIBRARIES})

```

The `main.cpp` file:

```cpp
#include <xfdtd/boundary/pml.h>
#include <xfdtd/monitor/field_monitor.h>
#include <xfdtd/monitor/movie_monitor.h>
#include <xfdtd/object/object.h>
#include <xfdtd/parallel/mpi_support.h>
#include <xfdtd/shape/sphere.h>
#include <xfdtd/simulation/simulation.h>
#include <xfdtd/waveform_source/tfsf_3d.h>

int main() {
  xfdtd::MpiSupport::setMpiParallelDim(
      1, 1, 2); // Set MPI parallel dim (1x1x2). If you don't want to use MPI,
                // it will do nothing.
  constexpr double dl {7.5e-3}; // set the grid resolution

  auto domain{std::make_shared<xfdtd::Object>(
      "domain",
      std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, -0.175, -0.175},
                                    xfdtd::Vector{0.35, 0.35, 0.35}),
      xfdtd::Material::createAir())}; // create calculation domain

  auto dielectric_sphere{std::make_shared<xfdtd::Object>(
      "dielectric_sphere",
      std::make_unique<xfdtd::Sphere>(xfdtd::Vector{0, 0, 0}, 0.1),
      std::make_unique<xfdtd::Material>(
          "a", xfdtd::ElectroMagneticProperty{3, 2, 0, 0}))};

  constexpr auto l_min{dl * 20};
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr std::size_t tfsf_start{
      static_cast<size_t>(15)}; // Set TFSF start index
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, 0, 0, 0,
      xfdtd::Waveform::gaussian(
          tau, t_0))}; // Add TFSF source with Gaussian waveform

  // Create movie monitor for Ex field in XZ plane
  auto movie_ex_xz{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(xfdtd::Vector{-0.175, 0, -0.175},
                                        xfdtd::Vector{0.35, dl, 0.35}),
          xfdtd::EMF::Field::EX, "", ""),
      10, "movie", "./data/dielectric_sphere_scatter/movie_ex_xz")};

  // Create simulation object with 2 threads
  auto s{xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{2, 1, 1}}};
  s.addObject(domain);
  s.addObject(dielectric_sphere);
  s.addWaveformSource(tfsf);
  // Add PML boundary
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  s.addMonitor(movie_ex_xz);
  s.run(2000); // Run simulation for 2000 steps
}

```

build the project:

```bash
cmake -B ./build
cmake --build ./build
```

run the project:

```bash
./build/xfdtd_first_example
```

If it goes well, you will see the output in the `./data/dielectric_sphere_scatter/movie_ex_xz` folder. As you had defined the `movie_ex_xz` monitor to save the Ex field in the XZ plane. The output will be a series of `.npy` files. You can use the `numpy` library to read and use the `matplotlib` library to plot the data in `Python`.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

if __name__ == '__main__':
    data_dir = './data/dielectric_sphere_scatter/movie_ex_xz'

    def get_sorted_files(dir: str):
        files = os.listdir(dir)
        files = [f for f in files if f.endswith('.npy')]
        files.sort()
        return [os.path.join(dir, f) for f in files]

    def read_data(file):
        return np.load(file)

    def process_data(data):
        d = np.squeeze(data)
        d[np.abs(d) < 1e-4] = 1e-4
        d = 10 * np.log10(np.abs(d))
        return d

    files = get_sorted_files(data_dir)

    data = read_data(files[0])
    print(data.shape)
    data = process_data(data)
    print(data.shape)

    f, ax = plt.subplots()
    im = ax.imshow(data, cmap='jet', vmin=-25, vmax=5)
    f.colorbar(im)

    def update(frame):
        im.set_data(process_data(read_data(files[frame])))
        ax.set_title(f'frame {frame}')
        return im,

    ani = animation.FuncAnimation(f, update, frames=len(files),
                        interval=100, repeat=True)
    
    ani.save('movie_ex_xz.gif', writer='ffmpeg', fps=10)

```

the gif file will be saved in the current directory.

![movie_ex_xz](./doc/image/movie_ex_xz.gif)

## Parallel Computing

There are mainly three levels of parallel computing: Vectorization, Shared Memory and Distributed Memory.
We use the C++ Standard Thread, MPI and CUDA for shared memory and use MPI for distributed memory.

### Use C++ Standard Thread

You can use the following code to set the thread number while creating the simulation object.

```cpp
#include <xfdtd/simulation/simulation.h>

// Create simulation object with 2 threads in X direction, 1 thread in Y direction and 1 thread in Z direction
auto s{xfdtd::Simulation{dl, dl, dl, 0.9, xfdtd::ThreadConfig{2, 1, 1}}};
```

We suggest that you set the thread dimension to 1 in the X and Y direction and set 2 or 4 in the Z direction.

### Use MPI

XFDTD CORE support MPI parallel computing. You can use the following command to compile the project with MPI.

First build XFDTD_CORE with MPI

In XFDTD_CORE project directory

```bash
cmake -DXFDTD_CORE_WITH_MPI=ON -B ./build
cmake --build ./build
cmake --build ./build --target install
```

And then you can return to your project directory and use the following `CMakeLists.txt` file to compile your project with MPI.

You need to link your target with `MPI::MPI_CXX` library.

```cmake
# You can check if the MPI is enabled by checking the XFDTD_CORE_WITH_MPI variable.
if(XFDTD_CORE_WITH_MPI)
  target_link_libraries(your_target PRIVATE MPI::MPI_CXX)
endif()
```

```bash
mpiexec -n ${num_of_core} ./build/your_executable
```

If you want to config the MPI parallel dim, you can use the following code.

```cpp
#include <xfdtd/parallel/mpi_support.h>

/*
some code
*/

// before you create the simulation
// set MPI parallel dim (2x2x1)
xfdtd::MpiSupport::setMpiParallelDim(2, 2, 1);
```

Then you can run the project with the following command.

```bash
mpiexec -n 4 ./build/your_executable
```

### CUDA

You can see the project [xfdtd_cuda](https://github.com/Mrwatermolen/XFDTD_CUDA) for the CUDA version of the XFDTD project.

## More Examples

See the `examples` folder for more examples.

[Example README](./examples/RAEDME.md)

## Known Issues

1. S11 parameter calculation is not correct.
2. Can't see the std::cout output in the console when using MPI in Linux.

## Acknowledgement

* Computational electrodynamics: the finite-difference time-domain method
* The Finite-Difference Time-Domain Method: Electromagnetics with MATLAB Simulations
* Parallel finite-difference time-domain method
* [uFDTD](https://github.com/john-b-schneider/uFDTD)
