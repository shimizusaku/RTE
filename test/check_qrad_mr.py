import numpy as np
import matplotlib.pyplot as plt



# p = np.load('/mnt/solar08b/c0287shimizu/work/data/d001/MultiRay/qrad_mr_test_python.npz')
# f_h = np.load('/mnt/solar08b/c0287shimizu/work/data/d001/MultiRay/qrad_mr_test_fortran.npz')
# f = np.load('/mnt/solar08b/c0287shimizu/work/data/d001/MultiRay/qrad_mr_test_fortran_notqrad.npz')
p = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/python.npz')
f = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/fortran.npz')

# get data
qrad_python = p['qrad']   # python Qrad data
# qrad_hotta = f_h['qrad_mr']
qrad_fortran = f['qrad']

# check difference
err_ = np.abs(qrad_fortran - qrad_python) / (qrad_fortran)
# err_ = np.abs(qrad_hotta - qrad_fortran) / (qrad_fortran + 1e-10)
# err_ = np.abs(qrad_hotta - qrad_python) / (qrad_hotta + 1e-10)
err = err_[:255, 1:256, 1:256]

ave_err = np.zeros(255)
for i in range(0, 255):
    ave_err[i] = np.average(err[i, :, :])
#     # print(ave_err[i])



# max err
max_err = np.amax(err)
max_err_index = np.where(err == max_err)
print('max_err      : ', max_err)
print('max_err_index: ', max_err_index)

# average err
# ave_err = np.average(err)
# print('average err: ', ave_err)


# plot
# plot = err[200, :, :]
plot = ave_err[:]
x = np.arange(255)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
# plt.imshow(plot)
plt.plot(plot, x)
# plt.title('qrad err')
plt.title('qrad average error')

# ticks = np.linspace(plot.min(), plot.max(), 10)
# cbar = plt.colorbar()
# cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/test/check'
file_name = 'qrad_err_mr.png'
file_path = os.path.join(save_dir, file_name)

# plt.show()
plt.savefig(file_path)
plt.close()