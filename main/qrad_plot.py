import numpy as np
import matplotlib.pyplot as plt

d = np.load('/mnt/solar07a/c0287shimizu/data/d001/MultiRay/fortran.npz')
qrad_tu1 = d['qrad_tu1']

plot = qrad_tu1[1:, 1:]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)
plt.title('multi directional')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)
cbar.set_label('Qrad')

save_dir = '/mnt/solar07a/c0287shimizu/fig/d001_fortran/MultiRay'
file_name = 'Qrad_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.tight_layout()

plt.savefig(file_path)
plt.close()