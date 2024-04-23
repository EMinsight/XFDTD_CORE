# XFDTD Core Examples

Some examples to demonstrate the usage of the XFDTD Core library.

## Getting Started

### Prerequisites

* C++ compiler supporting C++20
* CMake 3.20 or later
* XFDTD Core library

### Building

target list:

* cylinder_scatter_2d: 2D TEz scattering by a cylinder
* debye_sphere_scatter: 3D scattering by a Debye sphere
* drude_sphere_scatter: 3D scattering by a Drude sphere
* lorentz_sphere_scatter: 3D scattering by a Lorentz sphere
* dielectric_cube_scatter: 3D scattering by a dielectric cube
* dielectric_sphere_scatter: 3D scattering by a dielectric sphere
* half_wave_dipole: 3D half-wave dipole radiation
* inverted_f_antenna: 3D inverted-F antenna radiation
* microstrip_branch_line_coupler: 3D microstrip branch-line coupler
* microstrip_line: 3D microstrip line
* microstrip_low_pass: 3D microstrip low-pass filter
* quarter_wave_transformer: 3D quarter-wave transformer
* rlc_circuit: 3D RLC circuit
* simple_capacitor: 3D simple capacitor
* simple_circuit: 3D simple circuit
* slab_refection_transmission: 1D slab reflection and transmission

```bash
cmake -S . -B build
cmake --build build -t ${target}
```

If you want to build all targets, you can use the following command:

```bash
cmake --build build -t build-examples
```

Enable MPI support:

```bash
cmake -S . -B build -DXFDTD_CORE_WITH_MPI=ON
cmake --build build -t ${target}
```

### Running

all targets can be found in the `build/bin` directory.

```bash
./build/bin/${target}
```

Running with MPI:

```bash
mpirun -n ${n} ./build/bin/${target}
```

Running all targets:

```bash
chmod +x ./run_all_examples.sh
./run_all_examples.sh ./build/bin 
```

### Process Data

If nothing goes wrong, you can find the output data in the `./data` directory.

Use `Python` to process the data.

The script file in the `./plot_script` directory can be used to plot the data.

## Known Issues

Can't run with MPI.
