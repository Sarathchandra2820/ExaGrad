# ExaGrad

A Fortran implementation of Restricted Hartree-Fock (RHF) SCF for molecular electronic structure calculations, using [libcint](https://github.com/sunqm/libcint) for integral evaluation.

## Features

- RHF SCF with DIIS convergence acceleration
- Multiple Fock matrix construction methods:
  - `direct` — O(N^4) direct integral evaluation
  - `cholesky` — density-fitting with Cholesky decomposition
  - `block_cholesky` — blocked Cholesky density-fitting
  - `true_df` — true density-fitting with auxiliary basis
- OpenMP parallelization
- Electronic dipole moment calculation
- Automatic validation against PySCF reference energies

## Requirements

- gfortran
- LAPACK/BLAS (Accelerate on macOS, OpenBLAS on Linux)
- Python 3 with PySCF (for integral generation and validation)

## Build

```bash
make rhf_main
```

## Usage

```bash
./run.sh -g geometry/H2O.xyz -b cc-pVTZ -m true_df -a weigend -t 8
```

**Options:**

| Flag | Description | Example |
|------|-------------|---------|
| `-g` | Geometry file (XYZ) | `geometry/H2O.xyz` |
| `-b` | Basis set | `cc-pVTZ` |
| `-m` | Fock build method | `direct`, `cholesky`, `block_cholesky`, `true_df` |
| `-a` | Auxiliary basis (for DF methods) | `weigend` |
| `-t` | Number of OpenMP threads | `8` |

## Project Structure

```
src/
├── types/          # Molecule type definition
├── interfaces/     # libcint C bindings, BLAS/LAPACK wrappers
├── integrals/      # One-electron integral routines
├── scf/            # Fock builder and SCF driver
├── properties/     # Dipole moment integrals
└── programs/       # Main program entry point
geometry/           # Molecular geometries (XYZ format)
ints/               # Generated integral data files
python/             # PySCF integral export and validation scripts
```
