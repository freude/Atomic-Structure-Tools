import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


path='/home/mk/siesta_swarm/silicon_slab'

potential = Dataset(os.path.join(path, 'TotalPotential.grid.nc'))
charge = Dataset(os.path.join(path, 'TotalCharge.grid.nc'))
rho = Dataset(os.path.join(path, 'DeltaRho.grid.nc'))

print(charge.dimensions['xyz'])
a = ((np.array(potential.variables['gridfunc'][0, :, :, :])))
b = ((np.array(rho.variables['gridfunc'][0, :, :, :])))

size = a.shape
print(size)
plt.contour(a[0, :, :])
plt.show()
a1 = np.sum(a, axis=(0, 1)) / (size[0] * size[1])
b1 = np.sum(b, axis=(0, 1)) / (size[0] * size[1])
plt.plot(a[40, 40, :])
plt.plot(b[40, 40, :])
plt.show()
plt.plot(a1)
plt.plot(b1)
plt.show()

print(potential.variables['gridfunc'][0, :, :, 1])

from sisl import *
from sisl.io import fdfSileSiesta, bandsSileSiesta



from sisl.io.siesta import pdosSileSiesta

pdos = pdosSileSiesta(os.path.join(path, 'siesta.PDOS')).read_data()
geom = pdos[0]

from sisl import get_sile
fdf = get_sile(os.path.join(path, 'si_slab.fdf'))
geom = fdf.read_geometry()

energy = pdos[1]
pdos = pdos[2]

ldos = np.zeros((len(energy), len(geom)))
ind = 0

coords_min = [np.min(geom.xyz[:, 0]), np.min(geom.xyz[:, 1]), np.min(geom.xyz[:, 2])]
coords_max = [np.max(geom.xyz[:, 0]), np.max(geom.xyz[:, 1]), np.max(geom.xyz[:, 2])]


for j, g in enumerate(geom):
    for j1 in range(geom.orbitals[j]):
        ldos[:, j] += pdos[ind, :]
        ind += 1

args = geom[:, 0].argsort()
ldos = ldos[:, args]

from invdisttree import Invdisttree

xnew, ynew, znew = np.mgrid[-coords_min[0]:2*coords_max[0]:300j,
                            -coords_min[1]:2*coords_max[1]:10j,
                            -coords_min[2]:2*coords_max[2]:10j]

size = xnew.shape
ldos_grid = np.zeros((len(energy), size[0], size[1], size[2]))

for j in range(len(energy)):
    print(j)
    interp = Invdisttree(np.vstack((geom.xyz[:, 0], geom.xyz[:, 1], geom.xyz[:, 2])).T, ldos[j, :])
    ldos_grid[j, :, :, :] = interp(np.vstack((xnew.flatten(), ynew.flatten(), znew.flatten())).T,
                                   nnear=11,
                                   eps=0,
                                   p=1).reshape(xnew.shape)

print(ldos_grid)