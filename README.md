![q-e-logo](logo.jpg)

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO:
opEn-Source Package for Research in Electronic Structure, Simulation, and
Optimization)

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Repository description
This repository contains the experimental distribution of the Quantum ESPRESSO suite of codes
enabled for OpenMP offloading.

In [main](https://github.com/QEF/qe_test-omp5/tree/main) branch, the following QE components are
enabled for OpenMP Offload:

1.  FFTXlib
2. LAXlib
3. KS_Solvers
   * CG
   * Davidson
   * PPCG
   * ParO

plus some modules.

In [develop](https://github.com/QEF/qe_test-omp5/tree/devel) branch, upflib and PW app are also enabled for
OpenMP Offload.

## Usage
Both branches can be compiled with CMake only and need a new version of devXlib to be compiled separately.

In order to compile devXlib:

```
git clone https://gitlab.com/max-centre/components/devicexlib.git
cd devicexlib
git checkout devxlib_refactor
mkdir ./build
cd ./build
cmake -DCMAKE_PREFIX_PATH=$MKLROOT -DCMAKE_Fortran_COMPILER=ifx -DDEVXLIB_ENABLE_ACC=OPENMPGPU -DDEVXLIB_ENABLE_GPU_BLAS=MKLGPU -DCMAKE_INSTALL_PREFIX=<Choose an installation folder for the library> ..
make -j
make install
```

After devXlib has been compiled, this repository can be cloned and compiled using the following instructions:

```
git clone https://github.com/QEF/qe_test-omp5.git
cd qe_test-omp5
git checkout <main or devel>
mkdir ./build
cd ./build
cmake -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icc -DQE_ENABLE_OPENMP_OFFLOAD=yes -DDEVXLIB_ROOT=<devXlib installation folder that you chosen previously> ..
make -j
```

## Contributing
Quantum ESPRESSO is an open project: contributions are welcome.
Read the [Contribution Guidelines](CONTRIBUTING.md) to see how you
can contribute.

## LICENSE

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
