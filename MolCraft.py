from src.MoleculeData import POSCAR, AddMol
from src.AddMolecule import AddMolecule

# TODO: Add multiple molecules

# get input file path
input_file = input('Please input the POSCAR file name, defaul: [POSCAR]\n')
input_file = input_file if input_file.strip() else 'POSCAR'
input_addmol = input('Please input the AddMol file name, default: [AddMol]\n')
input_addmol = input_addmol if input_addmol.strip() else 'AddMol'

# check if the file exists
try:
    with open('./'+input_file, 'r'):
        pass
    with open('./'+input_addmol, 'r'):
        pass
except FileNotFoundError:
    print('File not found!')
    exit(1)

# read files
addmol = AddMol('./' + input_addmol)
poscar = POSCAR('./' + input_file)

def get_input_range(axis, default=[0, 1.0]):
    tmp = input(f'Please input the filled range of {axis}, default: [0.0 1.0]\n')
    if not tmp.strip():
        return default
    return [float(x) for x in tmp.split()]

# Get the xyz range of poscar in direct coordinates
x_range = get_input_range('x')
y_range = get_input_range('y')
z_range = get_input_range('z')
add_range = [x_range, y_range, z_range]

# check if the range is valid
try:
    assert 0 <= x_range[0] < x_range[1] <= 1
    assert 0 <= y_range[0] < y_range[1] <= 1
    assert 0 <= z_range[0] < z_range[1] <= 1
except AssertionError:
    print('Invalid range!')
    exit(1)

# Get the number of molecules you want to add
n_mol = input('Please input the number of molecules you want to add, default: [1]\n')
n_mol = int(n_mol) if n_mol.strip() else 1

# Get the cosntant for distance calculation
const_dist = input('Please input the constant for distance calculation, default: [0.5]\n')
const_dist = float(const_dist) if const_dist.strip() else 0.5

# Add molecules
if AddMolecule(poscar, add_range, n_mol, addmol, const_dist):
    print('Add molecules successfully!')
else:
    print('Add molecules failed!')

# Write the new POSCAR
poscar.to_direct()
poscar.write_POSCAR('./'+input_file+'_new')

