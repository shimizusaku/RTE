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


# average intensity
in_lc = np.zeros(255)
in_sc = np.zeros(255)

for i in range(255):
    in_lc[i] = np.average(intensity_lc_[i, :, :])
    in_sc[i] = np.average(intensity_sc[0, 1, 1, 1, i, :, :])
    # print(i, abs(in_lc[i] - in_sc[i]))

intensity_lc_m = np.zeros((255, 256, 256))
intensity_sc_m = np.zeros((255, 257, 257))

for i in range(255):
    intensity_lc_m[i, :, :] = intensity_lc_[i, :, :] - in_lc[i]
    intensity_sc_m[i, :, :] = intensity_sc[0, 1, 1, 1, i, :, :] - in_lc[i]
intensity_sc_m = intensity_sc_m[:, 1:, 1:]


# error
in_err_ = np.zeros((255, 256, 256))
in_err = np.zeros(255)

for i in range(255):
    for j in range(256):
        for k in range(256):
            in_err_[i, j, k] = abs(intensity_lc_m[i, j, k] - intensity_sc_m[i, j, k]) / intensity_sc_m[i, j, k]
            # print(in_err_[i, j, k])

for i in range(255):
    in_err[i] = np.average(in_err_[i, :, :])


# plot
x = np.arange(255)

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.plot(in_err[:], x)
plt.title('ave intensity error')
plt.xlabel('error')
plt.ylabel('height')

save_dir = '/mnt/solar07a/c0287shimizu/fig/test'
file_name = 'intensity_err.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()