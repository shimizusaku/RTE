import numpy as np
import matplotlib.pyplot as plt


# load
lc_data = np.load('/mnt/solar07a/c0287shimizu/data/d001/LC/test.npz')
intensity_lc = lc_data['intensity_p']
sc_data = np.load('/mnt/solar07a/c0287shimizu/data/d001/MultiRay/fortran.npz')
                  
intensity_sc = sc_data['intensity']


# average intensity
in_lc = np.zeros(255)
in_sc = np.zeros(255)

for i in range(255):
    in_lc[i] = np.average(intensity_lc[i, :, :])
    in_sc[i] = np.average(intensity_sc[0, 1, 1, 1, i, :, :])
    # print(i, abs(in_lc[i] - in_sc[i]))


# error
in_err = np.zeros(255)

for i in range(255):
    in_err[i] = abs(in_lc[i] - in_sc[i]) / in_sc[i]
    # print(i, in_err[i])


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