import numpy as np
import matplotlib.pyplot as plt


# load
d = np.load('/mnt/solar07a/c0287shimizu/data/d001/LC/test.npz')

i = d['intensity']

intensity = np.zeros(255)

for k in range(255):
    intensity[k] = np.average(i[0, :, :, k])

##### plot #####
print('# plot #')
plot = intensity[:]
x = np.arange(255)

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
# plt.imshow(plot)
plt.plot(plot, x)
plt.title('ave intensity lc')

save_dir = '/mnt/solar07a/c0287shimizu/fig/test'
file_name = 'lc.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()