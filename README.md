# MolCraft

## Description

MolCraft is a tool designed to add molecules to `POSCAR` files. It reads `POSCAR` and `AddMol` files, generating a new `POSCAR` file with the added molecules. MolCraft uses van der Waals radii to ensure legal distances between molecules.

Key features include:
- Auto-detects POSCAR (`.vasp`) and AddMol files in the current directory with interactive menus.
- Supports any crystal system, not just cubic.
- Supports choosing the range of positions to add molecules.
- Supports periodic boundary conditions.
- Supports multiple AddMol files and random structure generation.

MolCraft checks the distances between the added molecules and the atoms in the crystal structure. If the distance between the added molecules and the atoms in the crystal structure is less than the sum of the van der Waals radii times a constant (`const_dist`), MolCraft will not add the molecules.

More specifically, for a atom A and a atom B, only if the following condition is satisfied, MolCraft will add:

$$ Dis(A, B) <  (R_A + R_B) \times const\_dist $$

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

### Interactive Mode

MolCraft uses an interactive prompt. On startup, it automatically detects available files:

1. **POSCAR Selection** — A numbered menu of `.vasp` / `POSCAR*` files is shown. Select by number, or press Enter for the default.

2. **AddMol Selection** — Similarly, detected molecule files are listed. You can select one or multiple files (e.g. `1 3`), or type `all` to select all.

3. **Fill Ranges** — Enter the fractional coordinate ranges (0–1) for molecule placement on each axis. Press Enter to accept defaults.

4. **Molecule Count** — Number of molecules to add per AddMol file. If multiple AddMol files were selected, you can set different counts for each.

5. **Distance Constant** — Scaling factor for vdW radii overlap check (default: 0.5).

6. **Random Structures** — Number of independent random configurations to generate.

7. **Confirmation** — A summary is shown before execution. Type `n` to cancel, or Enter to proceed.

### Input Files

- **POSCAR**: The crystal structure file in VASP format (`.vasp` or `POSCAR*`).
- **AddMol**: The molecule definition file. Format below.

**AddMol file format:**
```
num_atoms
x y z element
```

Example (`AddMol`):
```
1
-4.428486593    0.907560162    0.010519803  O
```

For a multi-atom molecule (e.g. H2O):
```
3
x1 y1 z1 O
x2 y2 z2 H
x3 y3 z3 H
```

### Output

Generated structures are saved to `./output/{poscar_name}/` as individual files named `{poscar_name}_new_{N}`.

## Test

Test files are included in the root directory:

- `10h2o.vasp` — Empty crystal structure (cubic box, 15A).
- `AddMol` — Single oxygen atom molecule definition.
- `Add_H2O` — Water molecule definition.

Run `python MolCraft.py` and follow the interactive prompts. For example, to fill the box with 192 water molecules at default settings (adjust the fill range to taste):

```
Select POSCAR file:
  [1] 10h2o.vasp [default]
  [2] Ir_O.vasp
Choice (Enter=default #1):

Select AddMol file(s):
  [1] AddMol [default]
  [2] Add_H2O
  ...
Choice (Enter=default #1): 2

Number of molecules to add per AddMol [default: 1]: 192
...
Start? [Y/n]:
```

## License

This project is licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details.

## Contact Information

Xiang Liu — shawnbox202025@gmail.com
