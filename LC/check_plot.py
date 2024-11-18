import numpy as np
import matplotlib.pyplot as plt


# load
lc_data = np.load('/mnt/solar07a/c0287shimizu/data/d001/LC/test.npz')
intensity_lc = lc_data['intensity_p']
sc_data = np.load('/mnt/solar07a/c0287shimizu/data/d001/MultiRay/fortran.npz')   
intensity_sc = sc_data['intensity']

intensity_lc = intensity_lc[:, 1:, 1:]
intensity_lc_ = np.zeros((255, 256, 256))

for k in range(256):
    for j in range(256):

        if j < 96 and k >= 96:
            intensity_lc_[:, j, k] = intensity_lc[:, 159+j, k-96]
        elif j >= 96 and k < 96:
            intensity_lc_[:, j, k] = intensity_lc[:, j-96, 159+k]
        elif j < 96 and k < 96:
            intensity_lc_[:, j, k] = intensity_lc[:, 159+j, 159+k]
        else:
            intensity_lc_[:, j, k] = intensity_lc[:, j-96, k-96]




# plot
plot = intensity_lc_[254, :, :]
plot_ = intensity_sc[0, 1, 1, 1, 254, :, :]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot, vmin=plot.min(), vmax=plot.max())
plt.xlabel('y')
plt.ylabel('z')

ticks = np.linspace(plot_.min(), plot_.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)
cbar.set_label('intensity')

save_dir = '/mnt/solar07a/c0287shimizu/fig/test'
file_name = 'intensity_lc_.png'
file_path = os.path.join(save_dir, file_name)

plt.tight_layout()
plt.savefig(file_path)
plt.close()