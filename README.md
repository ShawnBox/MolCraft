# MolCraft

<!-- Version is defined in `src/__init__.py` (`__version__`) and `pyproject.toml`. -->

**Version: 1.0.0**

## Description

MolCraft is a Python CLI tool that adds molecules to `POSCAR` (VASP format) crystal structures. It reads `POSCAR` and `AddMol` files and generates a new `POSCAR` file with the added molecules, using van der Waals radii to enforce legal interatomic distances. The molecule placement uses random Monte Carlo sampling with periodic boundary conditions and supports any crystal system.

Key features include:

- Auto-detects POSCAR (`.vasp`) and AddMol files in the current directory with interactive menus.
- Supports any crystal system, not just cubic.
- Supports choosing the fractional-coordinate range where molecules are added.
- Supports periodic boundary conditions.
- Supports multiple AddMol files and random structure generation.

MolCraft checks the distance between each added molecule and the atoms already in the crystal structure. If the distance between an added molecule atom A and a crystal atom B is less than the sum of their van der Waals radii times a constant (`const_dist`), the placement is rejected. More precisely, a molecule is only added when every pair satisfies:

$$ Dis(A, B) \ge (R_A + R_B) \times const\_dist $$

where $Dis(A, B)$ is the distance between atom A and atom B, and $R_A$, $R_B$ are their van der Waals radii read from `VdW.ini` (unknown elements default to 2.0 Å).

## Requirements

- **Python ≥ 3.10** (the code uses PEP 585 / PEP 604 type hints such as `dict[str, float]` and `list | np.ndarray`).
- **numpy** (the only runtime dependency).

## Installation

```sh
git clone https://github.com/ShawnBox/MolCraft.git
cd MolCraft

# Option A — install the package (installs the `src` package + numpy):
pip install -e .

# Option B — install only the runtime dependency and run from the repo root:
pip install numpy
```

Run MolCraft from the repository root (the entry point is `MolCraft.py`):

```sh
python MolCraft.py
# print the version:
python MolCraft.py --version   # or -v
```

## Usage

### Interactive Mode

MolCraft uses an interactive prompt. On startup it automatically detects available files:

1. **POSCAR Selection** — A numbered menu of `.vasp` / `POSCAR*` files is shown. Select by number, or press Enter for the default.
2. **AddMol Selection** — Detected molecule files are listed. You can select one or multiple (e.g. `1 3`), or type `all` to select all.
3. **Fill Ranges** — Enter the fractional coordinate ranges (0–1) for molecule placement on each axis. Press Enter to accept defaults.
4. **Molecule Count** — Number of molecules to add per AddMol file. If multiple AddMol files were selected, you can set different counts for each.
5. **Distance Constant** — Scaling factor for the vdW radius overlap check (default: 0.5).
6. **Random Structures** — Number of independent random configurations to generate.
7. **Confirmation** — A summary is shown before execution. Type `n` to cancel, or Enter to proceed.

### Scripted / Non-interactive Mode

Run without the interactive prompts by supplying `--poscar` and `--addmol`:

```sh
python MolCraft.py --poscar Ir_O.vasp --addmol Add_O --count 5 --seed 42 --yes
```

| Option | Meaning |
| --- | --- |
| `--poscar FILE` | POSCAR file to use (required in scripted mode). |
| `--addmol FILE...` | One or more AddMol files (required in scripted mode). |
| `--count N...` | Molecules per AddMol file — a single value broadcast to all files, or one value per file. |
| `--range XLO XHI YLO YHI ZLO ZHI` | Fractional-coordinate fill range on each axis (default `0 1 0 1 0 1`). |
| `--const F` | vdW scaling constant (default `0.5`). |
| `--n-rand N` | Number of independent random configurations (default `1`). |
| `--seed N` | Random seed for reproducible placements (also applies to interactive mode). |
| `--out DIR` | Output directory (default `./output/<poscar>/`). |
| `--yes` | Skip the confirmation prompt. |
| `--version`, `-v` | Print the version and exit. |

Supplying `--poscar` or `--addmol` switches to scripted mode; running with no options starts the interactive prompts. Using a fixed `--seed` gives reproducible results.

### Input Files

- **POSCAR**: The crystal structure file in VASP format (`.vasp` or `POSCAR*`).
- **AddMol**: The molecule definition file. Format:

```
num_atoms
x y z element
```

Example (`AddMol`, a single oxygen atom):

```
1
-4.428486593    0.907560162    0.010519803  O
```

Example of a multi-atom molecule (e.g. H2O):

```
3
x1 y1 z1 O
x2 y2 z2 H
x3 y3 z3 H
```

### Output

Generated structures are saved to `./output/{poscar_name}/` as individual files named `{poscar_name}_new_{N}`.

## Test

Test files are included in the repository root:

- `10h2o.vasp` — Empty crystal structure (cubic box, 15 Å).
- `AddMol` — Single oxygen atom molecule definition.
- `Add_H2O` — Water molecule definition.

Run `python MolCraft.py` and follow the interactive prompts. For example, to fill the box with water molecules:

```text
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

## Versioning

MolCraft follows [Semantic Versioning](https://semver.org/). The version is defined in **two places that must stay in sync**:

- `src/__init__.py` — `__version__`
- `pyproject.toml` — `[project].version`

Add a `--version` / `-v` flag to the CLI (`python MolCraft.py --version`) to print the current version.

## License

This project is licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details.

## Contact Information

Xiang Liu — shawnbox202025@gmail.com
