'''
This script is used to modify box size of a POSCAR file. (in z direction)
'''


from src.MoleculeData import POSCAR, AddMol

poscar = POSCAR('./mp2723-4layers.vasp_new')

poscar.to_cartesian()

poscar.box[2][2] += 15.

poscar.to_direct()

poscar.write_POSCAR('./mp2723-Ni-random.vasp')