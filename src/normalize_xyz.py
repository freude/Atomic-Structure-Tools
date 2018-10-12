import os


def norm_xyz(path_to_file, cell, norm, save=False):

    with open(path_to_file, 'r') as data:
        data = data.read()

    data = iter(data.splitlines())
    out = []

    for line in data:
        line = line.strip()
        if line.startswith(('Si', 'H')):
            list_of_words = line.split()
            try:
                x_coord = float(list_of_words[1])
                y_coord = float(list_of_words[2])
                z_coord = float(list_of_words[3])
                line = '{} {}  {}  {}'.format(list_of_words[0],
                                              x_coord/(norm * cell[0]),
                                              y_coord/(norm * cell[1]),
                                              z_coord/(norm * cell[2]))
            except ValueError:
                pass

        out.append(line)

    out = "\n".join(out)

    if save:
        with open(path_to_file, "w") as out_file:
            out_file.write(out)

    return out


cif_starter = """#======================================================================

# CRYSTAL DATA

#----------------------------------------------------------------------

data_VESTA_phase_1


_chemical_name_common                  'New structure'
_cell_length_a                         5.453368
_cell_length_b                         5.453368
_cell_length_c                         5.453368
_cell_angle_alpha                      90
_cell_angle_beta                       90
_cell_angle_gamma                      90
_space_group_name_H-M_alt              'P 1'
_space_group_IT_number                 1

loop_
_space_group_symop_operation_xyz
   'x, y, z'

loop_
   _atom_site_label
   _atom_site_occupancy
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_adp_type
   _atom_site_B_iso_or_equiv
   _atom_site_type_symbol

"""


def xyz2cif(path_to_file, cell, lattice_const, save=True):

    elem = ('Si', 'H')

    # cett unit cell parameters

    data = iter(cif_starter.splitlines())
    out = []

    for line in data:
        line = line.strip()

        if line.startswith('_cell_length_a'):
            line = "_cell_length_a                         {}".format(cell[0] * lattice_const)

        if line.startswith('_cell_length_b'):
            line = "_cell_length_b                         {}".format(cell[1] * lattice_const)

        if line.startswith('_cell_length_c'):
            line = "_cell_length_c                         {}".format(cell[2] * lattice_const)

        out.append(line)

    # read xyz file

    with open(path_to_file, 'r') as data:
        data = data.read()

    data = iter(data.splitlines())

    for line in data:
        line = line.strip()
        if line.startswith(elem):
            list_of_words = line.split()

            try:
                x_coord = float(list_of_words[1])
                y_coord = float(list_of_words[2])
                z_coord = float(list_of_words[3])
                if list_of_words[0].startswith('H'):
                    line = 'H  1.0   {}    {}    {}    Biso 1.0   H'.format(x_coord,
                                                                            y_coord,
                                                                            z_coord)
                if list_of_words[0].startswith('Si'):
                    line = 'Si  1.0   {}    {}    {}    Biso 1.0   Si'.format(x_coord,
                                                                            y_coord,
                                                                            z_coord)
                out.append(line)

            except ValueError:
                pass

    out = "\n".join(out)

    if save:
        head, tail = os.path.split(path_to_file)
        name = tail.split('.')[0]
        with open(os.path.join(head, name + '.cif'), "w") as out_file:
            out_file.write(out)

    return out


def cif2fdf(filename):

    os.system(("~/Downloads/cif2cell-1.2.10/cif2cell "
               "{} "
               "--no-reduce "
               "-p Siesta").format(filename))


if __name__ == "__main__":

    # path = '/home/mk/TB_work/SiNW2.xyz'
    # cell = [10, 10, 1]
    #
    # out = norm_xyz(path, cell, 5.5, save=True)
    #
    # a_Si = 5.453368
    #
    # xyz2cif(path, cell, a_Si, save=True)
    cif2fdf('/home/mk/TB_work/SiNW2_round.cif')
