from src.MoleculeData import POSCAR, AddMol
from src.AddMolecule import AddMolecule
from src import __version__
import os
import glob
import sys
import random
import argparse


# ---- AddMol file detection (content-based, robust) ----

def _is_addmol_file(path: str) -> bool:
    """Return True if *path* looks like a valid AddMol file.

    An AddMol file's first line is the number of atoms (a positive integer),
    followed by that many lines of ``x y z element``.
    """
    try:
        with open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except (OSError, UnicodeDecodeError):
        return False

    if not lines:
        return False

    first = lines[0].strip()
    if not first.isdigit():
        return False
    n = int(first)
    if n < 1:
        return False

    coord_lines = [l.split() for l in lines[1:1 + n]]
    if len(coord_lines) < n:
        return False

    for parts in coord_lines:
        if len(parts) < 4:
            return False
        try:
            float(parts[0])
            float(parts[1])
            float(parts[2])
        except ValueError:
            return False
        if not parts[3]:
            return False
    return True


def find_addmol_files():
    """Find probable AddMol files in the current directory by content."""
    candidates = []
    for name in os.listdir('.'):
        path = os.path.join('.', name)
        if not os.path.isfile(path):
            continue
        if name.startswith('.'):
            continue
        if os.path.splitext(name)[1].lower() in {'.py', '.vasp'}:
            continue
        if name in {'.gitignore', 'LICENSE', 'README.md', 'CLAUDE.md', 'pyproject.toml'}:
            continue
        if _is_addmol_file(path):
            candidates.append(name)
    return sorted(candidates)


def find_poscar_files():
    """Find POSCAR and .vasp files in the current directory (deduplicated)."""
    return sorted(dict.fromkeys(glob.glob("POSCAR*") + glob.glob("*.vasp")))


# ---- Interactive helpers ----

def select_file(prompt, files, default_idx=0, multi=False):
    """Show a numbered menu for file selection. Supports single and multi-select."""
    if not files:
        return None

    print(f"\n{prompt}")
    for i, f in enumerate(files):
        marker = " [default]" if i == default_idx else ""
        print(f"  [{i+1}] {f}{marker}")
    if multi:
        print("  (e.g. '1 3' for multiple, 'all' for all)")

    while True:
        choice = input(f"Choice (Enter=default #{default_idx+1}): ").strip()

        if not choice:
            return files if multi else [files[default_idx]]

        if multi and choice.lower() == 'all':
            return files

        try:
            indices = [int(x) for x in choice.split()]
        except ValueError:
            print(f"  Please enter valid number(s) (1-{len(files)}).")
            continue

        selected = [files[i-1] for i in indices if 1 <= i <= len(files)]
        invalid = [i for i in indices if not (1 <= i <= len(files))]
        if invalid:
            print(f"  Invalid: {invalid}. Valid range: 1-{len(files)}")
        if selected:
            return selected if multi else selected[0]


def get_float(prompt, default):
    while True:
        val = input(f"{prompt} [default: {default}]: ").strip()
        if not val:
            return default
        try:
            return float(val)
        except ValueError:
            print("  Enter a valid number.")


def get_int(prompt, default):
    while True:
        val = input(f"{prompt} [default: {default}]: ").strip()
        if not val:
            return default
        try:
            return int(val)
        except ValueError:
            print("  Enter a valid integer.")


def get_range(axis, default=(0, 1.0)):
    lo, hi = default
    while True:
        val = input(f"Fill range on {axis}-axis [default: {lo} {hi}]: ").strip()
        if not val:
            return list(default)
        parts = val.split()
        if len(parts) != 2:
            print("  Enter two numbers (e.g. '0 1').")
            continue
        try:
            a, b = float(parts[0]), float(parts[1])
        except ValueError:
            print("  Enter two valid numbers.")
            continue
        if not (0 <= a < b <= 1):
            print("  Must satisfy 0 <= low < high <= 1.")
            continue
        return [a, b]


# ---- Shared execution ----

def _generate_structure(poscar, add_range, n_mol, addmol_objs, const_dist,
                        n_rand, base_name, output_dir):
    """Run the placement and write n_rand independent structures."""
    os.makedirs(output_dir, exist_ok=True)
    for i in range(n_rand):
        p = poscar.copy()
        p.to_cartesian()

        if AddMolecule(p, add_range, n_mol, addmol_objs, const_dist):
            print(f"  [OK] Structure #{i+1}: molecules added")
        else:
            print(f"  [FAIL] Structure #{i+1}: failed (max tries reached)")

        p.write_POSCAR(os.path.join(output_dir, f'{base_name}_new_{i+1}'))

    print(f"\nDone! Output -> {output_dir}")


# ---- Interactive mode ----

def run_interactive():
    print("=" * 60)
    print(f"  MolCraft {__version__} - Add molecules to POSCAR crystal structures")
    print("=" * 60)

    # -- File selection --
    poscar_files = find_poscar_files()
    addmol_files = find_addmol_files()

    # POSCAR
    if poscar_files:
        sel = select_file("Select POSCAR file:", poscar_files)
        input_file = sel if sel else input("\nPOSCAR file name: ").strip() or 'POSCAR'
    else:
        input_file = input("\nPOSCAR file name: ").strip()
    if not input_file or not os.path.isfile(input_file):
        print(f"Error: '{input_file}' not found. Exiting.")
        return
    print(f"  > POSCAR: {input_file}")

    # AddMol
    if addmol_files:
        sel = select_file("Select AddMol file(s):", addmol_files, multi=True)
        input_addmols = sel if sel else []
    else:
        input_addmols = []

    if not input_addmols:
        raw = input("AddMol file name(s) (space-separated): ").strip().split()
        input_addmols = raw if raw else ['AddMol']

    missing = [f for f in input_addmols if not os.path.isfile(f)]
    if missing:
        print(f"Error: files not found: {', '.join(missing)}. Exiting.")
        return
    print(f"  > AddMol: {', '.join(input_addmols)}")

    # -- Read data --
    addmol_objs = [AddMol(f'./{f}') for f in input_addmols]
    poscar = POSCAR(f'./{input_file}')

    # -- Ranges --
    print("\n--- Fill ranges (direct coordinates, 0-1) ---")
    x_range = get_range('x')
    y_range = get_range('y')
    z_range = get_range('z')
    add_range = [x_range, y_range, z_range]

    # -- Molecule counts --
    n_mol_default = get_int("\nNumber of molecules to add per AddMol", 1)
    if len(input_addmols) > 1:
        n_mol = [get_int(f"  > Molecules for '{f}'", n_mol_default)
                 for f in input_addmols]
    else:
        n_mol = [n_mol_default]

    # -- Other parameters --
    const_dist = get_float("Distance constant (vdW scaling)", 0.5)
    n_rand = get_int("Random structures to generate", 1)

    # -- Confirmation --
    print(f"\n{'=' * 60}")
    print(f"  POSCAR:           {input_file}")
    print(f"  AddMol:           {', '.join(input_addmols)}")
    print(f"  Fill range:       x=[{x_range[0]},{x_range[1]}]  "
          f"y=[{y_range[0]},{y_range[1]}]  z=[{z_range[0]},{z_range[1]}]")
    print(f"  Molecules/AddMol: {n_mol}")
    print(f"  vdW constant:     {const_dist}")
    print(f"  Structures:       {n_rand}")
    print(f"{'=' * 60}")

    if input("Start? [Y/n]: ").strip().lower() in ('n', 'no'):
        print("Cancelled.")
        return

    base_name = os.path.splitext(input_file)[0]
    output_dir = f'./output/{base_name}/'
    _generate_structure(poscar, add_range, n_mol, addmol_objs, const_dist,
                        n_rand, base_name, output_dir)


# ---- Scripted (non-interactive) mode ----

def build_parser():
    p = argparse.ArgumentParser(
        prog='molcraft',
        description='MolCraft — add molecules into POSCAR crystal structures.')
    p.add_argument('-v', '--version', action='version',
                   version=f'MolCraft {__version__}')
    p.add_argument('--seed', type=int, default=None,
                   help='Random seed for reproducibility (also applies to interactive mode).')
    p.add_argument('--poscar', help='POSCAR file to use (scripted mode).')
    p.add_argument('--addmol', nargs='+', help='AddMol file(s) to add (scripted mode).')
    p.add_argument('--range', nargs=6, type=float,
                   metavar=('XLO', 'XHI', 'YLO', 'YHI', 'ZLO', 'ZHI'),
                   default=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                   help='Fill range (direct coords) per axis, e.g. "0 1 0 1 0 1".')
    p.add_argument('--count', nargs='+', type=int,
                   help='Molecules per AddMol file (single value or one per file).')
    p.add_argument('--const', type=float, default=0.5,
                   help='vdW scaling constant (default 0.5).')
    p.add_argument('--n-rand', type=int, default=1,
                   help='Number of random configurations (default 1).')
    p.add_argument('--out', default=None,
                   help='Output directory (default ./output/<poscar>/).')
    p.add_argument('--yes', action='store_true',
                   help='Skip the confirmation prompt (scripted mode).')
    return p


def run_scripted(parser, args):
    if args.poscar is None or args.addmol is None:
        parser.error("'--poscar' and '--addmol' are both required in scripted mode.")
    if not os.path.isfile(args.poscar):
        parser.error(f"POSCAR file not found: {args.poscar}")
    missing = [f for f in args.addmol if not os.path.isfile(f)]
    if missing:
        parser.error(f"AddMol files not found: {', '.join(missing)}")

    add_range = [[args.range[0], args.range[1]],
                 [args.range[2], args.range[3]],
                 [args.range[4], args.range[5]]]
    for lo, hi in add_range:
        if not (0 <= lo < hi <= 1):
            parser.error("--range values must satisfy 0 <= low < high <= 1.")

    if args.count is None:
        n_mol = [1] * len(args.addmol)
    elif len(args.count) == 1:
        n_mol = [args.count[0]] * len(args.addmol)
    elif len(args.count) == len(args.addmol):
        n_mol = list(args.count)
    else:
        parser.error("--count must be a single value or one value per AddMol file.")

    n_rand = max(1, args.n_rand)

    addmol_objs = [AddMol(f'./{f}') for f in args.addmol]
    poscar = POSCAR(f'./{args.poscar}')

    print(f"\n{'=' * 60}")
    print(f"  POSCAR:           {args.poscar}")
    print(f"  AddMol:           {', '.join(args.addmol)}")
    print(f"  Fill range:       x=[{add_range[0][0]},{add_range[0][1]}]  "
          f"y=[{add_range[1][0]},{add_range[1][1]}]  "
          f"z=[{add_range[2][0]},{add_range[2][1]}]")
    print(f"  Molecules/AddMol: {n_mol}")
    print(f"  vdW constant:     {args.const}")
    print(f"  Structures:       {n_rand}")
    print(f"  Seed:             {args.seed}")
    print(f"{'=' * 60}")

    if not args.yes and input("\nStart? [Y/n]: ").strip().lower() in ('n', 'no'):
        print("Cancelled.")
        return

    base_name = os.path.splitext(args.poscar)[0]
    output_dir = args.out or f'./output/{base_name}/'
    _generate_structure(poscar, add_range, n_mol, addmol_objs, args.const,
                        n_rand, base_name, output_dir)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.seed is not None:
        random.seed(args.seed)

    if args.poscar is not None or args.addmol is not None:
        run_scripted(parser, args)
    else:
        run_interactive()


if __name__ == '__main__':
    main()
