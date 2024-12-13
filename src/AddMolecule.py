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
    add_count = 0

    for i in range(MAX_TRIES * num_mol):
        if add_count == num_mol:
            return True

        pos = np.dot(get_random_pos(add_range), poscar.box)
        angle = get_random_angle()
        coords = addmol.rotate(angle)
        atoms = []
        for j in range(addmol.num_mol):
            atoms.append(Atom(addmol.elements[j], coords[j] + pos, ['T', 'T', 'T']))

        # test code
        # for j in range(len(atoms)):
        #     print(atoms[j])

        if not poscar.add_molecule(atoms, const_dist):
            continue
        add_count += 1
        print(f'Successfully added {add_count} molecules when trying {i} times')

    return False
        



        
