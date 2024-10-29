import numpy as np
import matplotlib.pyplot as plt


# load
p = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/1Ray/python.npz')
f = np.load('/mnt/solar07a/c0287shimizu/work/data/d001/1Ray/fortran.npz')


# data
qrad_python  = p['qrad'] 
qrad_fortran = f['qrad']

# check difference
err_ = np.abs(qrad_fortran - qrad_python) / (qrad_fortran)
err = err_[:, 1:256, 1:256]

ave_err = np.zeros(255)
for i in range(0, 255):
    ave_err[i] = np.average(err[i, :, :])



# max err
max_err = np.amax(err)
max_err_index = np.where(err == max_err)
print('max_err      : ', max_err)
print('max_err_index: ', max_err_index)



# plot
plot = ave_err[:]
x = np.arange(255)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.plot(plot, x)
plt.title('qrad average error')


save_dir = '/mnt/solar07a/c0287shimizu/work/test/check'
file_name = 'qrad_err_1r.png'
file_path = os.path.join(save_dir, file_name)

# plt.show()
plt.savefig(file_path)
plt.close()