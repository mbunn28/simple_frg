import numpy as np
import h5py
from numpy import complex128, linalg

f = h5py.File('src/magnus2.h5', 'r')
U = f['U']
U = U[()]
V = f['V']
V = V[()]
hopping = f['hopping']
nk = f['nk']
nk = nk[()]
nkf = f['nkf']
nkf = nkf[()]
mu = hopping['mu']
mu = mu[()]
t = hopping['t']
t = t[()]
tp = hopping['tp']
tp = tp[()]

print(nkf)
print(nk)
print(np.shape(V))
print(mu)
print(t)
print(tp)

A = 1

tensor = np.zeros((2,2,2,2), dtype=complex128)
for o1 in range(2):
    for o2 in range(2):
        for o3 in range(2):
            for o4 in range(2):
                if (o1 == o3 and o2 == o4):
                    tensor[o1,o2,o3,o4] = A
                if (o1 == o4 and o2 == o3):
                    tensor[o1,o2,o3,o4] = -A

u, s, vh = linalg.svd(tensor)
print(u)
print(s)
print(vh)