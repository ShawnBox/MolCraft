from src.MoleculeData import POSCAR, AddMol
from src.AddMolecule import AddMolecule
import os

# TODO: Add multiple molecules

# get input file path
input_file = input('Please input the POSCAR file name, defaul: [POSCAR]\n')
input_file = input_file if input_file.strip() else 'POSCAR'
input_addmols = input('Please input the AddMol file name, default: [AddMol]\n')
input_addmols = input_addmols.split() if input_addmols.strip() else ['AddMol']

# check if the file exists
try:
    with open('./'+input_file, 'r'):
        pass
    for input_addmol in input_addmols:
        with open('./'+input_addmol, 'r'):
            pass
except FileNotFoundError:
    print('File not found!')
    exit(1)

# read files
addmol = []
for input_addmol in input_addmols:
    addmol.append(AddMol('./'+input_addmol))

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
n_mol = n_mol.split() if n_mol.strip() else ['1']
n_mol = [int(x) for x in n_mol]

# Get the cosntant for distance calculation
const_dist = input('Please input the constant for distance calculation, default: [0.5]\n')
const_dist = float(const_dist) if const_dist.strip() else 0.5

# 获取随机生成的结构数目
n_rand = input('Please input the number of random structures you want to generate, default: [1]\n')
n_rand = n_rand.split() if n_rand.strip() else ['1']

for i in range(int(n_rand[0])):
    p = poscar.copy()

    # Add molecules
    if AddMolecule(p, add_range, n_mol, addmol, const_dist):
        print('Add molecules successfully in the '+str(i)+'th structure!')
    else:
        print('Add molecules failed in the '+str(i)+'th structure!')

    # Write the new POSCAR
    p.to_direct()

    output_dir = './output/'+input_file+'/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    p.write_POSCAR(output_dir+input_file+'_new_'+str(i))


