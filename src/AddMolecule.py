import random
from src.MoleculeData import Atom, POSCAR, AddMol
import numpy as np
from typing import Sequence

MAX_TRIES = 1000


def get_random_pos(range_: Sequence[list[float]]) -> list[float]:
    return [random.uniform(range_[0][0], range_[0][1]),
            random.uniform(range_[1][0], range_[1][1]),
            random.uniform(range_[2][0], range_[2][1])]


def get_random_angle() -> list[float]:
    return [random.uniform(0, 360),
            random.uniform(0, 360),
            random.uniform(0, 360)]


def AddMolecule(poscar: POSCAR,
                add_range: Sequence[list[float]],
                num_mol: list[int],
                addmol: list[AddMol],
                const_dist: float) -> bool:
    """
    poscar: POSCAR object
    add_range: list of list, the range of x, y, z in direct coordinates
    num_mol: list of int, the number of molecules you want to add per AddMol
    addmol: list of AddMol object
    const_dist: float, the constant for distance calculation
    """
    for i in range(len(num_mol)):
        added = 0
        tries = 0
        single = addmol[i].single_atom

        while added < num_mol[i] and tries < MAX_TRIES:
            pos = np.dot(get_random_pos(add_range), poscar.box)
            angle = [0, 0, 0] if single else get_random_angle()
            coords = addmol[i].rotate(angle)

            atoms = [
                Atom(addmol[i].elements[k], coords[k] + pos, ['T', 'T', 'T'])
                for k in range(addmol[i].num_mol)
            ]

            if not poscar.add_molecule(atoms, const_dist):
                tries += 1
                continue

            added += 1
            tries += 1
            print(f'  [{added}/{num_mol[i]}] Added molecule from AddMol #{i+1}'
                  f' (trial {tries})')

        if added < num_mol[i]:
            return False

    return True
