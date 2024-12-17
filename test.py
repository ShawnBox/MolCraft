from src.MoleculeData import POSCAR, AddMol

poscar = POSCAR('./POSCAR_new')

poscar.to_cartesian()

poscar.box[2][2] += 15.

poscar.to_direct()

poscar.write_POSCAR('./POSCAR_10h2o')