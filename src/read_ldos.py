import numpy as np
from scipy.io import FortranFile

f1 = FortranFile('/home/mk/siesta_swarm/silicon_slab/si_siesta.LDOS1', 'r')


cell = f1.read_record(np.float)
sizes = f1.read_record(np.uint32)

ans = []

for j1 in range(sizes[1]*sizes[2]):
    ans.append(f1.read_record(np.float))

ans = np.array(ans)
ans = ans.reshape((108, 108, 320))
print('hi')





