import R2D2
import time
import rte_mr
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
d.read_qq(tstep, 'se')  

ro0 = d.p['ro0']                                                       # background field density
ro1 = d.qq['ro']                                                       # density disturbances
te0 = d.p['te0']                                                       # background field temperature
te1 = d.qq['te']                                                       # temperature disturbances
se0 = d.p['se0']                                                       # background field entropy
se1 = d.qq['se']                                                       # entropy disturbances
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
Te0, tmp, tmp = np.meshgrid(te0, d.p['y'], d.p['z'], indexing='ij')
Se0, tmp, tmp = np.meshgrid(se0, d.p['y'], d.p['z'], indexing='ij')

# input data
ro = Ro0 + ro1                                                         # density
te = Te0 + te1                                                         # temperature
se = Se0 + se1                                                         # entropy
op = d.qq['op']                                                        # opacity
x = d.p['x']                                                           # height of index i
y = d.p['y']                                                           # lenght of index j
z = d.p['z']                                                           # lenght of index k


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
intensity, tu_1ray, x_tu1, qrad, qrad_tu1 = rte_mr.rte_multiray(ro, te, se, op, x, y, z)

##### save data #####
print('# save data #')
np.savez('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/fortran.npz', \
        intensity=intensity, tu_1ray=tu_1ray, x_tu1=x_tu1, qrad=qrad, qrad_tu1=qrad_tu1)


##### plot #####
print('# plot #')

# intensity 
x_point = [1, lx-1]
mti = ['0', '1' , '2']
iti = ['0', '1']
jti = ['0', '1']
kti = ['0', '1']

for m in range(0, 3):
    for idir in range(0, 2):
        for jdir in range(0, 2):
            for kdir in range(0, 2):
                plot = intensity[m, idir, jdir, kdir, x_point[idir], :, :]

                plt.rcParams['xtick.direction'] = 'in'
                plt.rcParams['ytick.direction'] = 'in'

                plt.figure()
                plt.imshow(plot)
                plt.title(mti[m] + iti[idir] + jti[jdir] + kti[kdir])

                ticks = np.linspace(plot.min(), plot.max(), 10)
                cbar = plt.colorbar()
                cbar.set_ticks(ticks)

                save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_fortran/MultiRay'
                file_name = mti[m] + iti[idir] + jti[jdir] + kti[kdir] + '.png'
                file_path = os.path.join(save_dir, file_name)

                plt.savefig(file_path)
                plt.close()



# Qrad
plot = qrad_tu1[1:, 1:]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)
plt.title('Qrad tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_fortran/MultiRay'
file_name = 'Qrad_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()



end_time = time.time()
elapsed_time = end_time - start_time
print('Calculation time:', elapsed_time)