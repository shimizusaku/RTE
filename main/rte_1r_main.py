import R2D2
import time
import rte_1r
import numpy as np
import matplotlib.pyplot as plt

start_time = time.time()


##### load R2D2 data #####
tstep = 10
data_path = '/mnt/solar07a/c0287shimizu/R2D2/run/d001/data/'
d = R2D2.R2D2_data(data_path)
d.read_qq(tstep,'ro')    
d.read_qq(tstep,'op')     
d.read_qq(tstep, 'te')  

ro0 = d.p['ro0']                                                       # background field density
ro1 = d.qq['ro']                                                       # density disturbances
te1 = d.qq['te']                                                       # temperature disturbances
te0 = d.p['te0']                                                       # background of tenperature
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
Te0, tmp, tmp = np.meshgrid(te0, d.p['y'], d.p['z'], indexing='ij')

# input data
ro = Ro0 + ro1                                                         # density
te = Te0 + te1                                                         # temperature
op = d.qq['op']                                                        # opacity
x = d.p['x']                                                           # height of index i



#####  information #####
# MHD grid number
nx = 256
ny = 256
nz = 256
# RAD grid number
lx = nx - 1
ly = ny + 1
lz = nz + 1

##### RTE #####
intensity, tu_1ray, x_tu1, qrad, qrad_tu1 = rte_1r.rte_1ray(ro, te, op, x)

##### save data #####
print('# save data #')
np.savez('/mnt/solar07a/c0287shimizu/work/data/d001/1Ray/fortran.npz', \
        intensity=intensity, tu_1ray=tu_1ray, x_tu1=x_tu1, qrad=qrad, qrad_tu1=qrad_tu1)

##### plot #####
print('# plot #')

# intensity 
x_point = [1, lx-1]
iti = ['0', '1']


for idir in range(0, 2):

    plot = intensity[idir, x_point[idir], 1:, 1:]

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.figure()
    plt.imshow(plot)
    plt.title(iti[idir])

    ticks = np.linspace(plot.min(), plot.max(), 10)
    cbar = plt.colorbar()
    cbar.set_ticks(ticks)

    save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_fortran/1Ray'
    file_name = iti[idir] + '.png'
    file_path = os.path.join(save_dir, file_name)

    plt.savefig(file_path)
    plt.close()


# Qrad
plot = qrad_tu1[1:, 1:]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)
plt.title('Qrad 1ray tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_fortran/1Ray'
file_name = 'Qrad_1ray_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()




end_time = time.time()
elapsed_time = end_time - start_time
print('Calculation time:', elapsed_time)