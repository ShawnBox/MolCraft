import random
from src.MoleculeData import Atom
import numpy as np

# Constants
'''
Here, we define the maximum number of tries to add a molecule to the system.
Maybe we will change this value in the future.
'''
MAX_TRIES = 1000

def get_random_pos(range):
    x = random.uniform(range[0][0], range[0][1])
    y = random.uniform(range[1][0], range[1][1])
    z = random.uniform(range[2][0], range[2][1])
    return [x, y, z]

def get_random_angle():
    return [random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360)]

def AddMolecule(poscar, add_range, num_mol, addmol, const_dist):
    r"""
    poscar: POSCAR object
    add_range: list of list, the range of x, y, z in direct coordinates
    num_mol: list of int, the number of molecules you want to add
    addmol: list of AddMol object
    const_dist: float, the constant for distance calculation
    """

    try_count = 0
    for i in range(len(num_mol)):
        add_count = 0
        while True:
            if try_count == MAX_TRIES:
                return False
            if add_count == num_mol[i]:
                break
            pos = np.dot(get_random_pos(add_range), poscar.box)
            angle = get_random_angle()
            coords = addmol[i].rotate(angle)
            atoms = []
            for k in range(addmol[i].num_mol):
                atoms.append(Atom(addmol[i].elements[k], coords[k] + pos, ['T', 'T', 'T']))

            if not poscar.add_molecule(atoms, const_dist):
                continue
            add_count += 1
            try_count += 1
            print(f'Successfully added {add_count} molecules of the {i}th AddMol when trying {try_count} times')
        
    return True
        



        
