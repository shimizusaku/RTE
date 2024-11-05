import R2D2
import numpy as np
import test_lc
# import matplotlib.pyplot as plt
from math import log, exp, pi, copysign



##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('/mnt/solar07a/c0287shimizu/R2D2/run/d001/data/')
d.read_qq(10,'ro')    
d.read_qq(10,'op')    
d.read_qq(10, 'te')   

ro0 = d.p['ro0']      # background field density
ro1 = d.qq['ro']      # density disturbances
te1 = d.qq['te']       # temperature disturbances
te0 = d.p['te0']       # background of tenperature
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
Te0, tmp, tmp = np.meshgrid(te0, d.p['y'], d.p['z'], indexing='ij')
                                    

# input data
ro = Ro0 + ro1                                                         # density
te = Te0 + te1                                                         # temperature
op = d.qq['op']                                                        # opacity
x = d.p['x']                                                           # height of index i
y = d.p['y']                                                           # lenght of index j
z = d.p['z']                                                           # lenght of index k



##### RTE #####
intensity = test_lc.rte_longcharacteristic(ro, te, op, x, y, z)


##### save data #####
print('# save data #')
np.savez('/mnt/solar07a/c0287shimizu/data/d001/LC/test.npz', intensity=intensity)

# ##### plot #####
# print('# plot #')
# plot = tu_1ray[1, 100, 100, :]
# x = np.arange(257)

# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'

# plt.figure()
# # plt.imshow(plot)
# plt.plot(plot, x)
# plt.title('tu_1ray')

# save_dir = '/mnt/solar07a/c0287shimizu/work/LC/fig'
# file_name = 'tu_1ray.png'
# file_path = os.path.join(save_dir, file_name)

# plt.savefig(file_path)
# plt.close()