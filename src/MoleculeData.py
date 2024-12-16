import numpy as np

# read VdW.ini file
Vdw = {}
with open('VdW.ini', 'r') as f:
    lines = f.readlines()
    for line in lines:
        tmp = line.split()
        key = tmp[0]
        value = tmp[2]
        Vdw[key] = float(value)

class POSCAR():
    def __init__(self, file_path):
        self.file_path = file_path
        self.lattice_constant = 0
        self.box = np.zeros((3, 3))
        self.elements = []

        self.dynamic = False
        self.dynamic_label = []

        self.direct = True
        self.coordinates = []
        self.atoms = []
        self.element_dict = {}

        self.ex_atoms = []

        self.read_POSCAR()

    def read_POSCAR(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()
        
        # read lattice constant
        self.lattice_constant = float(lines[1].strip())

        # read box vectors
        for i in range(2, 5):
            self.box[i-2] = np.array([float(x) for x in lines[i].split()])

        # read elements
        tmp_name = lines[5].split()
        tmp_num = [int(lines[6].split()[k]) for k in range(len(tmp_name))]
        for i, element in enumerate(tmp_name):
            self.element_dict[element] = tmp_num[i]

        def get_elements_idx(elements):
            idx = []
            for key, value in elements.items():
                idx += [key] * value
            return idx
        self.elements = get_elements_idx(self.element_dict)

        # remove information has been read
        lines = lines[7:]

        if lines[0].strip():

            if  lines[0].split()[0].lower() in ['s', 'selective', 'dynamics', 'd']:
                self.dynamic = True
                lines = lines[1:]

            # read direct or cartesian
            if lines[0].strip().lower() in ['c', 'cartesian']:
                self.direct = False

            # read coordinates and dynamic label
            for i in range(1, 1 + len(self.elements)):
                self.coordinates.append([float(x) for x in lines[i].split()[:3]])
                if self.dynamic:
                    self.dynamic_label.append(lines[i].split()[3:])
            
        # make sure the coordinates are in cartesian
        self.to_cartesian()

        # get atoms
        for i in range(len(self.elements)):
            if self.dynamic:
                self.atoms.append(Atom(self.elements[i], self.coordinates[i], self.dynamic_label[i]))
            else:
                self.atoms.append(Atom(self.elements[i], self.coordinates[i]))

        # get extra atoms considering the periodic boundary conditions
        for atom in self.atoms:
            self.add_new_ex_atoms(atom)

    def to_cartesian(self):
        if not self.direct:
            return True
        self.direct = False
        for i in range(len(self.coordinates)):
            self.coordinates[i] = np.dot(self.coordinates[i], self.box)
        for atom in self.atoms:
            atom.coordinates = np.dot(atom.coordinates, self.box)
        for ex_atom in self.ex_atoms:
            ex_atom.coordinates = np.dot(ex_atom.coordinates, self.box)
        return True
    
    def to_direct(self):
        if self.direct:
            return True
        self.direct = True
        for i in range(len(self.coordinates)):
            self.coordinates[i] = np.dot(self.coordinates[i], np.linalg.inv(self.box))
        for atom in self.atoms:
            atom.coordinates = np.dot(atom.coordinates, np.linalg.inv(self.box))
        for ex_atom in self.ex_atoms:
            ex_atom.coordinates = np.dot(ex_atom.coordinates, np.linalg.inv(self.box))
        return True

    def add_molecule(self, new_atoms, const_dist=0.5):
        '''
        Here, we should note that when the number of atoms in the new molecule is large, the time complexity of this function is O(nm) (m is the number of ex_atom, n is the number of new_atoms), which is not efficient.
        Maybe we can use the k-d tree to optimize this function, so the time complexity can be reduced to O(nlogm).
        '''
        # TODO: using k-d tree to optimize this function
        for ex_atom in self.ex_atoms:
            for new_atom in new_atoms:
                if not ex_atom.is_legal(new_atom, const_dist):
                    return False

        self.atoms += new_atoms
        self.elements += [new_atom.element for new_atom in new_atoms]
        self.coordinates += [new_atom.coordinates for new_atom in new_atoms]
        
        for new_atom in new_atoms:
            self.add_new_ex_atoms(new_atom)

        if self.dynamic:
            self.dynamic_label += [new_atom.dynamic_label for new_atom in new_atoms]

        for new_atom in new_atoms:
            if new_atom.element in self.element_dict:
                self.element_dict[new_atom.element] += 1
            else:
                self.element_dict[new_atom.element] = 1
        return True

    def write_POSCAR(self, file_path):
        with open(file_path, 'w') as f:
            f.write('POSCAR generated by MolCraft\n')
            f.write(str(self.lattice_constant) + '\n')
            for i in range(3):
                f.write(f'\t{self.box[i][0]:.10f}\t{self.box[i][1]:.10f}\t{self.box[i][2]:.10f}\n')

            for key in self.element_dict.keys():
                f.write(f'{key} ')
            f.write('\n')
            for value in self.element_dict.values():
                f.write(f'{value} ')
            f.write('\n')

            if self.dynamic:
                f.write('Selective Dynamics\n')

            f.write('Direct\n')
            for key in self.element_dict.keys():
                for atom in self.atoms:
                    if atom.element == key:
                        f.write(f'\t{atom.coordinates[0]:.10f}\t{atom.coordinates[1]:.10f}\t{atom.coordinates[2]:.10f}')
                        if self.dynamic:
                            f.write(' ' + ' '.join(atom.dynamic_label))
                        f.write('\n')

    def add_new_ex_atoms(self, atom):
        def is_edge(coordinate):
            for i in range(3):
                if coordinate[i] < -0.3 or coordinate[i] > 1.3:
                    return False
            return True

        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    new_coordinates = atom.coordinates + np.dot(np.array([i, j, k]), self.box)
                    if not is_edge(np.dot(new_coordinates, np.linalg.inv(self.box))):
                        continue
                    self.ex_atoms.append(Atom(atom.element, new_coordinates))
            


class Atom():
    def __init__(self, element, coordinates, dynamic_label=[]):
        self.element = element
        self.coordinates = coordinates
        self.dynamic_label = dynamic_label
        self.Vdw = Vdw[element]
    
    def __str__(self):
        return f'{self.element}: {self.coordinates[0]} {self.coordinates[1]} {self.coordinates[2]} {self.dynamic_label}'

    def is_legal(self, atom, const_dist=0.5):
        if np.linalg.norm(np.array(self.coordinates) - np.array(atom.coordinates)) < (self.Vdw + atom.Vdw) * const_dist:
            return False
        return True


class AddMol():
    def __init__(self, file_path):
        self.file_path = file_path
        self.num_mol = 0
        self.elements = []
        self.coordinates = []

        self.read_AddMol()
        self.geometric_center = self.geometric_center()
        self.coordinates = self.new_coordinates(self.geometric_center)
    
    def read_AddMol(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()
        
        self.num_mol = int(lines[0].strip())

        for i in range(1, 1 + self.num_mol):
            self.coordinates.append([float(x) for x in lines[i].split()[:3]])
            self.elements.append(lines[i].split()[3])

    # to calculate the geometric center of the added molecules
    def geometric_center(self):
        return np.mean(self.coordinates, axis=0)
    
    # to calculate new coordinates of the added molecules
    def new_coordinates(self, center):
        return np.array(self.coordinates) - center
    
    # to calculate new coordinates after rotation
    def rotate(self, angle_list):
        angle_list = np.array(angle_list) * np.pi / 180
        matrix_x = np.array([[1, 0, 0], [0, np.cos(angle_list[0]), -np.sin(angle_list[0])], [0, np.sin(angle_list[0]), np.cos(angle_list[0])]])
        matrix_y = np.array([[np.cos(angle_list[1]), 0, np.sin(angle_list[1])], [0, 1, 0], [-np.sin(angle_list[1]), 0, np.cos(angle_list[1])]])
        matrix_z = np.array([[np.cos(angle_list[2]), -np.sin(angle_list[2]), 0], [np.sin(angle_list[2]), np.cos(angle_list[2]), 0], [0, 0, 1]])
        matrix = np.dot(matrix_x, np.dot(matrix_y, matrix_z))
        return np.dot(self.coordinates, matrix)


if __name__ == '__main__':
    poscar = POSCAR('./POSCAR')
    poscar.to_cartesian()
    poscar.write_POSCAR('./POSCAR_cartesian')
    poscar.to_direct()
    poscar.write_POSCAR('./POSCAR_direct')
