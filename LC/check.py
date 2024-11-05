import numpy as np
import matplotlib.pyplot as plt


# load
d = np.load('/mnt/solar07a/c0287shimizu/data/d001/LC/test.npz')

intensity = d['intensity']

##### plot #####
print('# plot #')
plot = intensity[1, 100, 100, ]
x = np.arange(256)

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
# plt.imshow(plot)
plt.plot(plot, x)
plt.title('tu_1ray_sc')

save_dir = '/mnt/solar07a/c0287shimizu/work/LC/fig'
file_name = 'sc.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()