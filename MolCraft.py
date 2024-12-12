from src.MoleculeData import POSCAR, AddMol
from src.AddMolecule import AddMolecule

# read files
addmol = AddMol('./AddMol')
poscar = POSCAR('./POSCAR')

def get_input_range(axis):
    tmp = input(f'Please input the filled range of {axis}, like: 0 1.0\n')
    return [float(x) for x in tmp.split()]

# Get the xyz range of poscar in direct coordinates
x_range = get_input_range('x')
y_range = get_input_range('y')
z_range = get_input_range('z')
add_range = [x_range, y_range, z_range]

# Get the number of molecules you want to add
n_mol = int(input('Please input the number of molecules you want to add\n'))

# Add molecules
if AddMolecule(poscar, add_range, n_mol, addmol):
    print('Add molecules successfully!')
else:
    print('Add molecules failed!')

# Write the new POSCAR
poscar.to_direct()
poscar.write_POSCAR('./POSCAR_new')

