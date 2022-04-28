import numpy as np
from numpy import ndarray
import h5py
from numpy import complex128, linalg
import matplotlib.pyplot as plt

file_name = "adaeuler"
f = h5py.File(f'src/{file_name}.h5', 'r')
# print(f.keys())
U = f['U']
U = U[()]
reV = f['reV']
reV = reV[()]
imV = f['imV']
imV = imV[()]
V = reV + 1j*imV
V = np.roll(V, (1,1,1,1,1,1), (0,1,2,3,4,5))
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

LamC = 0.05

print("Superconducting Gap Calculation")
print("--------------------------------")
print(f"nkf={nkf} (not used)")
print(f"nk={nk}")
print(f"t={t}")
print(f"tp={tp}")
print(f"mu={mu}")
print(f"lambda_crit={LamC}")


def ferm(x):
    return 1/(1+np.exp(x))

def e(kx,ky):
    return -2*t*(np.cos(kx)+np.cos(ky))-4*tp*np.cos(kx)*np.cos(ky) - mu

def L(kpx, kpy):
    return (ferm(-e(kpx,kpy)/LamC) - ferm(e(kpx,kpy)/LamC))/(2 * e(kpx,kpy))
    
VL = np.zeros((nk*nk,nk*nk), dtype=complex128)

for kx in range(nk):
    for ky in range(nk):
        for kpx in range(nk):
            for kpy in range(nk):
                VL[kx+nk*ky, kpx+nk*kpy] = V[kpx,kpy,kx,ky,-kx,-ky] * L(kpx,kpy)

u, s, vh = linalg.svd(VL)

del_vec = u[:,0]
# del_vec = vh[0,:]
delta = del_vec.reshape((nk,nk))
# print(delta)

plt.pcolormesh(np.real(delta),cmap="coolwarm")
plt.colorbar()
plt.show()
#plt.savefig(f"scgap_{file_name}.png")

# Gamma = np.zeros((nk,nk),dtype=complex128)
# for i in range(nk):
#     for j in range(nk):
#         Gamma[i,j] = V[0,0,i,j,-i,-j] - V[i,j,0,0,0,0]

# plt.clf()
# plt.imshow(np.abs(Gamma))
# plt.show()
