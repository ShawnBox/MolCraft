from src.MoleculeData import POSCAR, AddMol
from src.AddMolecule import AddMolecule
import os
import glob


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


def find_addmol_files():
    """Find probable AddMol files in the current directory."""
    all_items = os.listdir('.')
    excluded_ext = {'.py', '.vasp', '.md', '.txt'}
    excluded_names = {'.gitignore', 'LICENSE', 'README.md'}
    candidates = []
    for name in all_items:
        path = os.path.join('.', name)
        if not os.path.isfile(path):
            continue
        ext = os.path.splitext(name)[1].lower()
        if ext in excluded_ext or name.startswith('.') or name in excluded_names:
            continue
        candidates.append(name)
    return sorted(candidates)


def find_poscar_files():
    """Find POSCAR and .vasp files in the current directory."""
    return sorted(glob.glob("POSCAR*") + glob.glob("*.vasp"))


def main():
    print("=" * 60)
    print("  MolCraft - Add molecules to POSCAR crystal structures")
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
        n_mol = []
        for i, f in enumerate(input_addmols):
            n_mol.append(get_int(f"  > Molecules for '{f}'", n_mol_default))
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

    # -- Run --
    base_name = os.path.splitext(input_file)[0]
    output_dir = f'./output/{base_name}/'
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


if __name__ == '__main__':
    main()
