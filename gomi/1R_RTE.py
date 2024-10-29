import os
import time
import R2D2
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, pi

start_time = time.time()

##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('/mnt/solar07a/c0287shimizu/R2D2/run/d002/data/')
d.read_qq(10,'ro')   
d.read_qq(10,'op')   
d.read_qq(10, 'te')  

ro0 = d.p['ro0']                   # background field density
ro1 = d.qq['ro']                   # density disturbances
op = d.qq['op']                    # optical depth
x = d.p['x']                       # height of index i
y = d.p['y']                       # lenght of index j
z = d.p['z']                       # lenght of index k
te1 = d.qq['te']                    # temperature disturbances
te0 = d.p['te0']                    # background of tenperature



#####  information #####
# MHD grid number
nx = 254
ny = 254
nz = 254
# RAD grid number
lx = nx - 1
ly = ny + 1
lz = nz + 1


# ray direction
step = [-1, 1]                  # step direction
shift = [-1, 0]                 # shift of index
shift_s = [0, -1]               # not shift 
ista = [lx-1, 1]                # start point of i loop
iend = [1, lx-1]                # finish point of i loop
mux_1ray = 1 /sqrt(3)



# list
alpha = np.zeros((lx, ly, lz))                                 # absorption_coefficient in RAD grid
beta = np.zeros((lx, ly, lz))                                    # planck function in RAD grid
delta_tu = np.zeros((2, lx, ly, lz))                             # delta optical depth
exp_delta_tu = np.zeros((2, lx, ly, lz))
tu_1ray = np.zeros((lx, ly, lz))                                   # optical depth
x_tu1 = np.zeros((ly, lz))
tu0 = 0.1    
delta_I = np.zeros((2, lx, ly, lz))                              # delta intensity
I = np.zeros((2, lx, ly, lz))                                    # intensity
Qrad = np.zeros((lx, ly, lz))                                 # radiative heating rate
Qrad_tu1 = np.zeros((ly, lz))



##### define #####
# absorption cefficient for RAD grid
def absorption_coefficient_RAD(alpha_MHD, i, j, k):
    return exp((log(alpha_MHD[i,   j,   k]) + log(alpha_MHD[i,   j,   k+1]) \
              + log(alpha_MHD[i,   j+1, k]) + log(alpha_MHD[i,   j+1, k+1]) \
              + log(alpha_MHD[i+1, j,   k]) + log(alpha_MHD[i+1, j,   k+1]) \
              + log(alpha_MHD[i+1, j+1, k]) + log(alpha_MHD[i+1, j+1, k+1])) \
            / 8)

# planck function for RAD grid
def planck_function_RAD(beta_MHD, i, j, k):
    return exp((log(beta_MHD[i,   j,   k]) + log(beta_MHD[i,   j,   k+1]) \
              + log(beta_MHD[i,   j+1, k]) + log(beta_MHD[i,   j+1, k+1]) \
              + log(beta_MHD[i+1, j,   k]) + log(beta_MHD[i+1, j,   k+1]) \
              + log(beta_MHD[i+1, j+1, k]) + log(beta_MHD[i+1, j+1, k+1])) \
            / 8)

# delta optical depth
def rte_delta_tu(alpha_down, alpha_up, lray):
    # parameter to check
    epc = 1e-4

    alpha_du = alpha_down / alpha_up - 1
    # for zero divide
    alpha_dum = copysign(1, alpha_du) * max(abs(alpha_du), epc) 

    alpha_du15 = 1 - 0.5 * alpha_du
    # for zero divide
    alpha_du15m = copysign(1, alpha_du15) * max(abs(alpha_du15), epc) 

    delta_tu0 = alpha_up * alpha_dum / log(alpha_dum + 1)
    delta_tu1 = alpha_up / alpha_du15m
    
    tmp = 0.5 + copysign(0.5, abs(alpha_du) - epc)
    delta_tu = (delta_tu0 * tmp + delta_tu1 * (1 - tmp)) * lray

    return delta_tu

# delta intensity
def rte_delta_I(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu):
    # parameter to check
    epc = 1e-4

    exp_delta_tu = rte_exp(delta_tu)

    beta_upt = beta_up * exp_delta_tu

    beta_du = beta_upt / beta_down - 1
    beta_dum = copysign(1, beta_du) * max(abs(beta_du), epc) 

    beta_du15 = 1 - 0.5 * beta_du
    # for zero divide
    beta_du15m = copysign(1, beta_du15) * max(abs(beta_du15), epc) 

    tmp = 0.5 + copysign(0.5, abs(beta_du) - epc)
    log_beta_du_delta_tu0 = log_beta_down - log_beta_up + delta_tu
    log_beta_du_delta_tu1 = 1                                       # value is not important

    log_beta_du_delta_tu = \
        log_beta_du_delta_tu0 * tmp + log_beta_du_delta_tu1 * (1 - tmp)

    delta_in0 = (beta_down - beta_upt) / log_beta_du_delta_tu
    delta_in1 = beta_down / beta_du15m

    delta_in = (delta_in0 * tmp + delta_in1 * (1 - tmp)) * delta_tu
    
    return delta_in

# def rte_exp(delta_tu, idir, i, j, k):
#     xmin = -700
#     tmp = 0.5 + copysign(0.5, - delta_tu[idir, i, j, k] - xmin)

#     return exp(max(- delta_tu[idir, i, j, k], xmin)) * tmp
def rte_exp(delta_tu):
    xmin = -700
    tmp = 0.5 + copysign(0.5, - delta_tu - xmin)

    return exp(max(- delta_tu, xmin)) * tmp






##### absorption coefficient & planck function #####
print('##### alpha & planck function #####')

Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
ro = Ro0 + ro1                                                        # density

alpha_data = ro * op                                                   # alpha_MHD_data is no Periodic Boundary Conditions 
alpha_MHD_ = np.insert(alpha_data, 0, alpha_data[:, ny-1, :], axis = 1)
alpha_MHD_ = np.insert(alpha_MHD_, ny+1, alpha_data[:, 0, :], axis = 1)
alpha_MHD = np.insert(alpha_MHD_, 0, alpha_MHD_[:, :, nz-1], axis = 2)   # alpha in MHD grid (Periodic Boundary Conditions)
alpha_MHD = np.insert(alpha_MHD, nz+1, alpha_MHD_[:, :, 0], axis = 2)


Te0 , tmp, tmp = np.meshgrid(te0, d.p['y'], d.p['z'], indexing = 'ij')
te = Te0 + te1                                                           # tempurature
sig = 5.67 * 10**(-5)                                                 # Stefan-Boltzmann constant

beta_data = (sig * te**4) / pi                                            # beta_MHD_data is no Periodic Boundary Conditions 
beta_MHD_ = np.insert(beta_data, 0, beta_data[:, ny-1, :], axis = 1)
beta_MHD_ = np.insert(beta_MHD_, ny+1, beta_data[:, 0, :], axis = 1)
beta_MHD = np.insert(beta_MHD_, 0, beta_MHD_[:, :, nz-1], axis = 2)            # B in MHD grid (Periodic Boundary Conditions)
beta_MHD = np.insert(beta_MHD, nz+1, beta_MHD_[:, :, 0], axis = 2)


for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            alpha[i, j, k] = absorption_coefficient_RAD(alpha_MHD, i, j, k)
            beta[i, j, k] = planck_function_RAD(beta_MHD, i, j, k)





##### intensity #####
print('##### intensity #####')

for idir in range(2):
    print('idir: ', idir)

    for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
        for j in range(0, ly):
            for k in range(0, lz):
                
                # downstream
                alpha_down     = alpha[i + shift[idir], j, k]
                log_alpha_down = log(alpha_down)

                beta_down     = beta[i + shift[idir], j, k]
                log_beta_down = log(beta_down)


                # upstream
                alpha_up = alpha[i + shift_s[idir], j, k]
                log_alpha_up = log(alpha_up)

                beta_up = beta[i + shift_s[idir], j, k]
                log_beta_up = log(beta_up)


                # ray lenght lray(m) in short characteristics
                lray = (x[i] - x[i-1]) / mux_1ray 

                
                # delta optical depth & delta intensity
                delta_tu[idir, i, j, k] = rte_delta_tu(alpha_down, alpha_up, lray)
                # exp_delta_tu[idir, i, j, k] = rte_exp(delta_tu[idir, i, j, k], idir, i, j, k)
                exp_delta_tu[idir, i, j, k] = rte_exp(delta_tu[idir, i, j, k])

                delta_I[idir, i, j, k] = rte_delta_I(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu[idir, i, j, k])

    
    # 1ray optical depth
    tu_1ray[lx-1, :, :] = 0
    x_tu1 = lx-1

    for i in range(lx-1, -1, -1):
        for j in range(0, ly):
            for k in range(0, lz):
                tu_1ray[i-1, j, k] = tu_1ray[i, j, k] + delta_tu[1, i, j, k]

                if (1 - tu_1ray[i, j, k]) * (1 - tu_1ray[i-1, j, k]) <= 0:
                    x_tu1[j, k] = i
    

    # intensity
    for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
        for j in range(0, ly):
            for k in range(0, lz):
                I[idir, i+shift[idir], j, k] = I[idir, i+shift_s[idir], j, k] * exp_delta_tu[idir, i, j, k] + delta_I[idir, i, j, k]





##### radiative heating rate #####
print('##### radiative heating rate #####')

for i in range(1, lx):
    for j in range(1, ly):
        for k in range(1, lz):
            # qrad_f
            fx_up   = (I[1, i,   j, k] - I[0, i,   j, k]) * mux_1ray * pi * 2
            fx_dpwn = (I[1, i-1, j, k] - I[0, i-1, j, k]) * mux_1ray * pi * 2

            # qrad
            Qrad[i, j, k] = - (fx_up - fx_dpwn) / ((x[i] - x[i-1]) * ro[i, j-1, k-1] * te[i, j-1, k-1])


for j in range(1, ly):
    for k in range(1, lz):
        Qrad_tu1[j, k] = Qrad[x_tu1[j, k], j, k]
                



##### save data #####
print('### save data ###')
np.savez('/mnt/solar07a/c0287shimizu/work/data/d001/1Ray/python.npz', qrad=Qrad, qrad_tu1=Qrad_tu1, intensity=I, tu1=tu, x_tu1=x_tu1)




##### plot #####
print('### plot ###')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# intensity
x_point = [1, lx-1]
iti = ['0', '1']

for i in range(0, 2):
    plot = I[x_point[i], :, :]

    plt.figure()
    plt.imshow(plot)
    plt.title(iti[i])

    ticks = np.linspace(plot.min(), plot.max(), 10)
    cbar = plt.colorbar()
    cbar.set_ticks(ticks)

    save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_py/1Ray'
    file_name = iti[i] + '.png'
    file_path = os.path.join(save_dir, file_name)

    plt.savefig(file_path)
    plt.close
    

# Qrad
plot = qrad_tu1[1:, 1:]

plt.figure()
plt.imshow(plot)
plt.title('Qrad_1ray tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_py/1Ray'
file_name = 'Qrad_1ray_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close





end_time = time.time()
elapsed_time = end_time - start_time
print('Calculation time:', elapsed_time)