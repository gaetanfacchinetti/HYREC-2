# HYREC-2

[![Build Status](https://github.com/gaetanfacchinetti/HYREC-2/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/gaetanfacchinetti/HYREC-2/actions/workflows/python-package.yml?query=branch%3Amaster)

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyhyrec)](https://pypi.org/project/pyhyrec/)
[![PyPi Version](https://img.shields.io/pypi/v/pyhyrec.svg)](https://pypi.python.org/pypi/pyhyrec/)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![Static Badge](https://img.shields.io/badge/physics-cosmology-darkblue)](https://en.wikipedia.org/wiki/Cosmology)
[![Static Badge](https://img.shields.io/badge/physics-recombination-yellow)](https://en.wikipedia.org/wiki/Recombination_(cosmology))



A code for primordial hydrogen and helium recombination including radiative transfer. An updated version of HyRec with new computing mode.
HYREC-2 is implemented in `CLASS` code as a default module for recombination history since version 3.0.0.
HYREC-2 can be used in `CAMB` code as well.

The original code HyRec developed by Yacine Ali-Haimoud and Chris Hirata (2010-17).
HYREC-2 developed with contributions from Nanoom Lee (2020).

Version: July 2020.
See readme.pdf for detailed info, and check HYREC-2 paper (Reference 1 below) for the code development.

## Using the python wrapper

In this fork we have implement a lightweight python wrapper for HYREC-2 using cython.

### Installation

* The easiest option for installation is to use
  ```
  pip install hyrec
  ```

### How to run HYREC-2

An example of notebook using the python wrapper is given [here](https://github.com/gaetanfacchinetti/HYREC-2/blob/master/example.ipynb).




## Using the classical C version

### Installation

* Clone the repo
   ```sh
   git clone https://github.com/nanoomlee/HYREC-2.git
   ```
### How to run HYREC-2

* Compiling the code
  ```sh
  gcc -lm -O3 hyrectools.c helium.c hydrogen.c history.c energy_injection.c hyrec.c -o hyrec
  ```
* Computing a recombination history
  ```sh
  ./hyrec < input.dat
  ```
The code will take parameters from `input.dat` and compute the recombination history. It will print out three columns in `output.dat`, which are: redshift, free electron fraction, matter temperature in Kelvin. For input parameter description, see `input.dat`.

## References

If you use HYREC-2, please cite the companion papers:

1. [N. Lee & Y. Ali-Haimoud, Phys. Rev. D 102 8, 083517 (2020)](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.102.083517)
2. [Y. Ali-Haimoud & C. M. Hirata, Phys. Rev. D 83, 043513 (2011)](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.83.043513)

You may also refer to the following papers, on which this work relies substantially:

1. [Y. Ali-Haimoud & C. M. Hirata, Phys. Rev. D 82, 063521 (2010)](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.82.063521)
2. [C. M. Hirata, Phys. Rev. D 78, 023001 (2008)](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.78.023001)
3. [E. R. Switzer & C. M. Hirata, Phys. Rev. D 77, 083008 (2008)](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.77.083008)

