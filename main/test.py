import numpy as np


#####  information #####
# MHD grid number
nx = 256
ny = 256
nz = 256
# RAD grid number
lx = nx - 1
ly = ny + 1
lz = nz + 1


p = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/test_python.npz')
f = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/test.npz')

Qrad_tu1 = np.zeros((ly, lz))

# tu_1ray = p['tu1_mr']
# qrad = p['qrad_mr']
tu_1ray = f['tu1_mr']
qrad = f['qrad_mr']

for j in range(1, ly):
    for k in range(1, lz):
        for i in range(lx-1, 0, -1):
            if (1 - tu_1ray[i, j, k])*(1 - tu_1ray[i-1, j, k]) <= 0:
                Qrad_tu1[j, k] = qrad[i, j, k]


# Qrad
plot = Qrad_tu1[1:, 1:]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)
plt.title('Qrad tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/test/fig'
file_name = 'Qrad_tu1_fortran.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()