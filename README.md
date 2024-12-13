# MolCraft

## Description

MolCraft is a tool designed to add molecules to `POSCAR` files. It reads `POSCAR` and `AddMol` files, generating a new `POSCAR` file with the added molecules. MolCraft uses van der Waals radii to ensure legal distances between molecules.

Key features include:
- Supports any crystal system, not just cubic.
- Supports choosing the range of positions to add molecules.
- Supports periodic boundary conditions.

MolCraft is suitable for researchers and engineers who need to add molecules to crystal structures.

MolCraft checks the distances between the added molecules and the atoms in the crystal structure. If the distance between the added molecules and the atoms in the crystal structure is less than the sum of the van der Waals radii times a constant (`const_dist`), MolCraft will not add the molecules.

More specifically, for a atom A and a atom B, only if the following condition is satisfied, MolCraft will add:

$$ Dis(A, B) <  (R_A + R_B) * const\underline{}dist $$

, where $Dis(A, B)$ is the distance between atom A and atom B, $R_A$ and $R_B$ are the van der Waals radii of atom A and atom B read from `VdW.ini`.

## Installation

Before installing MolCraft, make sure you have `Python 3.6` or later installed.

To install and set up MolCraft, follow these steps:

1. **Clone the repository:**
   ```sh
   git clone https://github.com/ShawnBox/MolCraft.git
   cd MolCraft
   ```
2. **Create a virtual environment (optional but recommended):**

   You can use `conda` or `venv` to create a virtual environment. 

3. **Install the required packages:**
   ```sh
   pip install numpy
   ```

4. **Run MolCraft:**
   ```sh
   python MolCraft.py
   ```

## Usage

1. **Input Files:**

   MolCraft requires two input files: `POSCAR` and `AddMol`.

   - `POSCAR`: The crystal structure file in VASP format.
   - `AddMol`: The molecule file in XYZ format.

   Please notice that if you want to add molecules to **an empty crystal structure**, you need to create a `POSCAR` file with the lattice parameters. And make sure that the `POSCAR` file has 9 lines even if line 6 to line 9 are empty in this case.

2. **User Input:**

   MolCraft will ask for the following information:

   - The file name of the `POSCAR` file.
   - The file name of the `AddMol` file.
   - The ranges of positions to add molecules.
   - The number of molecules to add.
   - The constant for distance calculation (`const_dist`).

## Test

Here, I add `POSCAR_H2O` and `AddMol_H2O` in the root directory as test input. `POSCAR_H2O` is a empty crystal structure file in VASP format, and `AddMol_H2O` is a water molecule file in XYZ format.

You can just run `python MolCraft.py` and input `POSCAR_H2O` and `AddMol_H2O`. Then you can use the default values for the ranges of positions to add molecules, and $192$ for the number of molecules to add, and $0.5$ for the constant for distance calculation.

It sinces that this version of MolCraft runs successfully but slowly. At least it was calculated in a minute and the output is:

```bash
Successfully added 192 molecules when trying 313 times
```

Maybe I will optimize the code in the future.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Contact Information

Xiang Liu - shawnbox202025@gmail.com
