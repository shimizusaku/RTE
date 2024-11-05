import numpy as np




d = np.load('/mnt/solar07a/c0287shimizu/data/d001/MultiRay/fortran.npz')

intensity = d['intensity']

inten = np.zeros(255)

for k in range(255):
    inten[k] = np.average(intensity[0, 1, 1, 1, k, :, :])

##### plot #####
print('# plot #')
plot = inten[:]
x = np.arange(255)

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
# plt.imshow(plot)
plt.plot(plot, x)
plt.title('ave intensity sc')

save_dir = '/mnt/solar07a/c0287shimizu/fig/test'
file_name = 'sc.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()