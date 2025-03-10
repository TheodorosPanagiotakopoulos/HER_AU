# Material Science Library for Hydrogen Evolution Reaction (HER) Simulations

## Overview
The tppy Library is designed to facilitate material science simulations related to the Hydrogen Evolution Reaction (HER). This library is built to work seamlessly with the Atomic Simulation Environment (ASE) and VASP. It provides a range of functions to manipulate atomic structures, extract properties, and prepare input files for simulations.

## Features
- **Atomic Manipulations:** Add, delete, swap, and shift atoms within structures.
- **Collective Coordinate Setup:** Create ICONST files for collective coordinates.
- **Charge Density Calculations:** Read and analyze CHGCAR files from VASP.
- **Structure Rotation and Movement:** Rotate atomic groups, move molecules symmetrically, and adjust atomic positions.
- **File Processing:** Extract useful information from VASP output files (OUTCAR), including Fermi energy and electronic charge.
- **Molecular Customization:** Add atoms symmetrically and adjust atomic bonds.

## Requirements
- **Python Version:** 3.7 or higher
- **Dependencies:**
  - `numpy`
  - `ase`

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/theodorosP/HER_AU.git
