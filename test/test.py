import numpy as np
import matplotlib.pyplot as plt


# load
f = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/test_python.npz')

qrad_tu1 = f['qrad_tu1_mr']
plot = qrad_tu1[1:, 1:]

plt.figure()
plt.imshow(plot)
plt.title('Qrad_1ray tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_py/MultiRay'
file_name = 'Qrad_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()