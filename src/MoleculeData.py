import numpy as np
import copy
import os
from typing import Optional

# ---- VdW radii loading ----

_Vdw: Optional[dict[str, float]] = None


def _load_vdw() -> dict[str, float]:
    """Lazy-load van der Waals radii from config file."""
    global _Vdw
    if _Vdw is not None:
        return _Vdw
    _Vdw = {}
    path = os.path.join(os.path.dirname(__file__), 'config', 'VdW.ini')
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    _Vdw[parts[0]] = float(parts[2])
    except (FileNotFoundError, ValueError):
        pass  # proceed with empty dict
    return _Vdw


def _get_vdw(element: str) -> float:
    return _load_vdw().get(element, 2.0)


# ---- Helper: PBC-aware minimum image distance ----

def _minimum_image_diff(coords_a: np.ndarray, coords_b: np.ndarray,
                        box: np.ndarray, inv_box: np.ndarray) -> np.ndarray:
    """
    Vectorised minimum-image distance vectors.
    coords_a: (N, 3) or (3,)   coords_b: (M, 3) or (3,)
    Returns (N, M, 3) or (N or M, 3) distance vectors.
    """
    a, b = np.atleast_2d(coords_a), np.atleast_2d(coords_b)
    diff_cart = a[:, None] - b[None, :]          # (N, M, 3)
    diff_frac = diff_cart @ inv_box               # (N, M, 3)
    diff_frac -= np.round(diff_frac)              # wrap to [-0.5, 0.5]
    return diff_frac @ box                        # (N, M, 3)


# ---- Classes ----

class Atom:
    def __init__(self, element: str, coordinates: list | np.ndarray,
                 dynamic_label: list | None = None):
        self.element = element
        self.coordinates = np.asarray(coordinates, dtype=float)
        self.dynamic_label = list(dynamic_label or [])
        self.Vdw = _get_vdw(element)

    def __str__(self) -> str:
        return (f'{self.element}: {self.coordinates[0]} '
                f'{self.coordinates[1]} {self.coordinates[2]} '
                f'{self.dynamic_label}')

    def is_legal(self, atom: 'Atom', const_dist: float = 0.5,
                 box: np.ndarray | None = None,
                 inv_box: np.ndarray | None = None) -> bool:
        """Check whether *atom* is too close (respecting PBC if *box* given)."""
        if box is not None and inv_box is not None:
            diff = _minimum_image_diff(self.coordinates, atom.coordinates,
                                       box, inv_box)
            dist = np.linalg.norm(diff.ravel())
        else:
            dist = np.linalg.norm(self.coordinates - atom.coordinates)
        return dist >= (self.Vdw + atom.Vdw) * const_dist


class POSCAR:
    """A VASP POSCAR crystal structure.

    ``self.atoms`` is the single source of truth for coordinates and element
    composition. ``coordinates`` and ``elements`` are derived properties, so
    the two representations never drift apart. Coordinates are kept in the
    representation indicated by ``self.direct``; call :meth:`to_cartesian` /
    :meth:`to_direct` before operating on them.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.lattice_constant = 0.0
        self.box = np.zeros((3, 3))
        self.dynamic = False
        self.direct = True
        self.atoms: list[Atom] = []
        self.element_dict: dict[str, int] = {}
        self._inv_box: np.ndarray | None = None

        self.read_POSCAR()

    # ---- Derived views (single source of truth is self.atoms) ----

    @property
    def coordinates(self) -> list[np.ndarray]:
        return [a.coordinates for a in self.atoms]

    @property
    def elements(self) -> list[str]:
        return [a.element for a in self.atoms]

    # ---- I/O ----

    def read_POSCAR(self) -> None:
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.lattice_constant = float(lines[1].strip())

        for i in range(3):
            self.box[i] = np.array([float(x) for x in lines[i + 2].split()])
        self._inv_box = np.linalg.inv(self.box)

        # elements
        names = lines[5].split()
        counts = [int(lines[6].split()[k]) for k in range(len(names))]
        self.element_dict = dict(zip(names, counts))
        elements = [e for e, n in self.element_dict.items() for _ in range(n)]

        lines = lines[7:]
        if not lines or not lines[0].strip():
            return

        if lines[0].split()[0].lower() in ('s', 'selective', 'dynamics', 'd'):
            self.dynamic = True
            lines = lines[1:]

        if lines and lines[0].strip().lower() in ('c', 'cartesian'):
            self.direct = False

        coords_raw, labels_raw = [], []
        for i in range(1, 1 + len(elements)):
            if i >= len(lines):
                break
            parts = lines[i].split()
            coords_raw.append(np.array([float(x) for x in parts[:3]]))
            labels_raw.append(parts[3:] if self.dynamic else [])

        self.atoms = [
            Atom(el, coord, label)
            for el, coord, label in zip(elements, coords_raw, labels_raw)
        ]

        self.to_cartesian()  # normalise to Cartesian on load

    def write_POSCAR(self, file_path: str) -> None:
        self.to_direct()
        with open(file_path, 'w') as f:
            f.write('POSCAR generated by MolCraft\n')
            f.write(str(self.lattice_constant) + '\n')
            for i in range(3):
                f.write(f'\t{self.box[i][0]:.10f}\t{self.box[i][1]:.10f}\t'
                        f'{self.box[i][2]:.10f}\n')
            for key in self.element_dict:
                f.write(f'{key} ')
            f.write('\n')
            for val in self.element_dict.values():
                f.write(f'{val} ')
            f.write('\n')
            if self.dynamic:
                f.write('Selective Dynamics\n')
            f.write('Direct\n')
            for key in self.element_dict:
                for atom in self.atoms:
                    if atom.element == key:
                        f.write(f'\t{atom.coordinates[0]:.10f}'
                                f'\t{atom.coordinates[1]:.10f}'
                                f'\t{atom.coordinates[2]:.10f}')
                        if self.dynamic:
                            f.write(' ' + ' '.join(atom.dynamic_label))
                        f.write('\n')

    # ---- Coordinate transforms ----

    def to_cartesian(self) -> bool:
        if not self.direct:
            return False
        self.direct = False
        for atom in self.atoms:
            atom.coordinates = atom.coordinates @ self.box
        return True

    def to_direct(self) -> bool:
        if self.direct:
            return False
        self.direct = True
        inv = self._inv_box if self._inv_box is not None else np.linalg.inv(self.box)
        for atom in self.atoms:
            atom.coordinates = atom.coordinates @ inv
        return True

    # ---- Core logic ----

    def add_molecule(self, new_atoms: list[Atom], const_dist: float = 0.5) -> bool:
        """
        Add *new_atoms* if none are within the vdW-radius distance threshold
        of existing atoms (respecting periodic boundary conditions).
        Uses vectorised numpy to check all pairs at once.
        """
        self.to_cartesian()

        if not self.atoms:
            self._append_atoms(new_atoms)
            return True

        existing = np.array([a.coordinates for a in self.atoms])   # (N, 3)
        existing_vdw = np.array([a.Vdw for a in self.atoms])       # (N,)
        inv = self._inv_box if self._inv_box is not None else np.linalg.inv(self.box)

        for new_atom in new_atoms:
            diff = _minimum_image_diff(existing, new_atom.coordinates,
                                       self.box, inv).ravel()       # (N,)
            dists = np.linalg.norm(diff.reshape(-1, 3), axis=1)    # (N,)
            limits = (existing_vdw + new_atom.Vdw) * const_dist
            if np.any(dists < limits):
                return False

        self._append_atoms(new_atoms)
        return True

    def _append_atoms(self, new_atoms: list[Atom]) -> None:
        self.atoms += new_atoms
        for a in new_atoms:
            self.element_dict[a.element] = self.element_dict.get(a.element, 0) + 1

    def copy(self) -> 'POSCAR':
        return copy.deepcopy(self)


class AddMol:
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.num_mol = 0
        self.elements: list[str] = []
        self.coordinates: np.ndarray = np.empty((0, 3))
        self.single_atom = False

        self.read_AddMol()
        centre = self.coordinates.mean(axis=0)
        self.coordinates = self.coordinates - centre

    def read_AddMol(self) -> None:
        with open(self.file_path, 'r') as f:
            lines = f.readlines()
        self.num_mol = int(lines[0].strip())
        self.single_atom = (self.num_mol == 1)
        coords = []
        for i in range(1, 1 + self.num_mol):
            parts = lines[i].split()
            coords.append([float(x) for x in parts[:3]])
            self.elements.append(parts[3])
        self.coordinates = np.array(coords, dtype=float)

    def rotate(self, angle_list: list[float]) -> np.ndarray:
        """Rotate coordinates by Euler angles (degrees) around geometric centre."""
        a = np.radians(angle_list)
        cx, sx = np.cos(a[0]), np.sin(a[0])
        cy, sy = np.cos(a[1]), np.sin(a[1])
        cz, sz = np.cos(a[2]), np.sin(a[2])
        Rx = np.array([[1, 0, 0],
                       [0, cx, -sx],
                       [0, sx,  cx]])
        Ry = np.array([[cy, 0, sy],
                       [0, 1, 0],
                       [-sy, 0, cy]])
        Rz = np.array([[cz, -sz, 0],
                       [sz,  cz, 0],
                       [0, 0, 1]])
        return self.coordinates @ (Rx @ Ry @ Rz)
