from __future__ import print_function, division
import numpy as np
import scipy.spatial
import matplotlib.pyplot as plt


def distance_to_surface(coords, surface):

    return np.inner(coords - surface, surface) / np.linalg.norm(surface)


def belong_to_surface(coords, surface):

    dist = distance_to_surface(coords, surface)
    flag = False

    if np.abs(dist) < 0.01:
        flag = True

    return flag


class KdTree(object):

    def __init__(self, coords, nn_dist=2.39):
        self._kd_tree = scipy.spatial.cKDTree(coords, leafsize=10)
        self._nn_distance = nn_dist

    def get_neighbours(self, coord):

        ans = self._kd_tree.query(coord,
                                  k=5,
                                  distance_upper_bound=self._nn_distance)

        ans1 = [ans[1][0]]

        for item in zip(ans[0], ans[1]):
            if self._nn_distance * 0.01 < item[0] < self._nn_distance:
                ans1.append(item[1])

        return ans1[1:]


def R_x(degrees):

    theta = np.pi / 180 * degrees

    R = np.matrix([[1.0, 0.0, 0],
                   [0.0, np.cos(theta), -np.sin(theta)],
                   [0.0, np.sin(theta), np.cos(theta)]])

    return R


def R_matrix(axis, degrees):

    theta = np.pi / 180 * degrees

    R = np.matrix(np.zeros((3, 3)))

    R[0, 0] = (1 - np.cos(theta)) * axis[0] ** 2 + np.cos(theta)
    R[0, 1] = (1 - np.cos(theta)) * axis[0] * axis[1] - axis[2] * np.sin(theta)
    R[0, 2] = (1 - np.cos(theta)) * axis[0] * axis[2] + axis[1] * np.sin(theta)
    R[1, 0] = (1 - np.cos(theta)) * axis[1] * axis[0] + axis[2] * np.sin(theta)
    R[1, 1] = (1 - np.cos(theta)) * axis[1] ** 2 + np.cos(theta)
    R[1, 2] = (1 - np.cos(theta)) * axis[1] * axis[2] - axis[0] * np.sin(theta)
    R[2, 0] = (1 - np.cos(theta)) * axis[2] * axis[0] - axis[1] * np.sin(theta)
    R[2, 1] = (1 - np.cos(theta)) * axis[2] * axis[1] + axis[0] * np.sin(theta)
    R[2, 2] = (1 - np.cos(theta)) * axis[2] ** 2 + np.cos(theta)

    return R


def add_atoms(coord, coords, length=0.7):

    a = []

    for item in coords:

        if len(coords) == 2:
            a1 = np.matrix(length * (coord - item))
            a1 = R_x(90) * a1.T
            a.append(coord + np.squeeze(np.array(a1)))

        if len(coords) == 1:
            a0 = np.matrix(length * (coord - item))
            a1 = R_x(90) * a0.T
            a1 = coord + np.squeeze(np.array(a1))
            a.append(a1)

            a1 = R_x(180) * R_x(90) * a0.T
            a1 = coord + np.squeeze(np.array(a1))
            a.append(a1)

    return a


def make_canted(coord, coords, atoms):

    atoms_out = []
    origin = np.array([coords[0][0], coord[1], coord[2]])
    axis = np.abs(coords[0] - coord)
    axis[0] = 0.0
    axis /= np.linalg.norm(axis)
    axis = np.squeeze(np.array(R_x(90) * np.matrix(axis).T))

    coord_out = origin + np.squeeze(np.array(R_matrix(axis, 10.0) * np.matrix(coord - origin).T))

    for atom in atoms:
        a1 = origin + np.squeeze(np.array(R_matrix(axis, 10.0) * np.matrix(atom - origin).T))
        atoms_out.append(a1)

    return coord_out, atoms_out


def fold_into(a1, cell):

    for j in range(3):
        if a1[j] < 0:
            a1[j] += cell[j]
        if a1[j] > cell[j]:
            a1[j] -= cell[j]

    return a1


def passivate_surface_cif(filename, elem, plane):

    coords = []
    atom_list = []
    coords_cart = []
    cell = np.array([0, 0, 0], dtype=float)

    flag = False

    try:
        data = open(filename, 'r')
        data_str = data.read()
    except IOError:
        data_str = filename

    header_str = data_str.split("_atom_site_type_symbol\n", 1)[0] + "_atom_site_type_symbol\n"

    data = iter(data_str.splitlines())

    for line in data:
        line = line.strip()

        if line.startswith('_cell_length_a'):
            cell[0] = float(line.split()[1])

        if line.startswith('_cell_length_b'):
            cell[1] = float(line.split()[1])

        if line.startswith('_cell_length_c'):
            cell[2] = float(line.split()[1])

        if line.startswith(elem):
            coord = np.array(line.split()[2:5], dtype=float)
            coords.append(coord)
            coords_cart.append(coord * cell)
            atom_list.append(line + "\n")

    coords = np.array(coords)

    for j in range(3):
        plane[j] = plane[j](coords[:, j])

    coords_cart = np.array(coords_cart)

    kd_tree = KdTree(coords_cart)

    for j, coord in enumerate(coords_cart):
        if belong_to_surface(coords[j], plane):
            ans = kd_tree.get_neighbours(coord)

            atoms = add_atoms(coord, coords_cart[ans])
            coord, atoms = make_canted(coord, coords_cart[ans], atoms)

            # for j in range(len(atoms)):
            #     atoms[j] = fold_into(atoms[j], cell)
            #
            # coord = fold_into(coord, cell)

            atom_list[j] = 'Si   1.0  {}   {}   {}  Biso  1.0   Si \n'.format(coord[0]/cell[0],
                                                                              coord[1]/cell[1],
                                                                              coord[2]/cell[2])

            for atom in atoms:
                a = atom.tolist()
                atom_list.append('H   1.0  {}   {}   {}  Biso  1.0   H \n'.format(a[0]/cell[0],
                                                                                  a[1]/cell[1],
                                                                                  a[2]/cell[2]))

    for atom in atom_list:
        header_str += atom

    with open('test.cif', 'w') as file:
        file.write(header_str)

    return header_str

#
# def fdf2cif(path):
#
#     with open(path, 'r') as data:


if __name__ == '__main__':

    filename_si_bulk = '/home/mk/Downloads/AMS_DATA.cif'
    filename_si_supercell = 'test.cif'

    import os
    os.system(("~/Downloads/cif2cell-1.2.10/cif2cell "
               "{} "
               "--supercell-translation-vector=[0.125,0,0] "
               "--supercell-postvacuum-translation=[0.05,0,0] "
               "--no-reduce "
               "--supercell=[3,1,5] "
               "--supercell-vacuum=[1,0,0] "
               "-p CASTEP").format(filename_si_bulk))

    def zero(x):
        return np.min(x * 0)

    cif = passivate_surface_cif(filename_si_supercell, ('Si', 'H'), [np.max, zero, zero])
    passivate_surface_cif(cif, ('Si', 'H'), [np.min, zero, zero])

    os.system(("~/Downloads/cif2cell-1.2.10/cif2cell "
               "{} "
               "--no-reduce "
               "-p Siesta").format(filename_si_supercell))
