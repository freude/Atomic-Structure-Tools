import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


path = '/home/mk/siesta_swarm/silicon_nw2'

data = Dataset(os.path.join(path, 'siesta.MD.nc'))

# for j in range()

coords = data.variables['xa']

num_steps, num_atoms, _ = data.variables['xa'].shape

lines = []

atom_dict = {'1': 'H', '2': 'Si'}

for j in range(num_steps):

    lines.append('{}\n'.format(num_atoms))

    for jj in range(num_atoms):

        x = str(data.variables['xa'][j, jj, 0].data)
        y = str(data.variables['xa'][j, jj, 1].data)
        z = str(data.variables['xa'][j, jj, 2].data)

        string = '{}    {}    {}    {}'.format(atom_dict[str(data.variables['isa'][jj])], x, y, z)
        lines.append(string)


lines = "\n".join(lines)

with open(os.path.join(path, 'frames.xyz'), "w") as out_file:
    out_file.write(lines)


import MDAnalysis


a = MDAnalysis.coordinates.XYZ.XYZReader(os.path.join(path, 'frames.xyz'))


print('hi')