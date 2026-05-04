'''
This script is used to modify box size of a POSCAR file. (in z direction)
'''


from src.MoleculeData import POSCAR, AddMol

poscar = POSCAR('./POSCAR')

poscar.to_cartesian()

poscar.box[2][2] -= 5.

poscar.to_direct()

poscar.write_POSCAR('./POSCAR_new')