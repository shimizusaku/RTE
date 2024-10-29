# This code is RTE for multi rays.
import os
import time
import R2D2
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, pi

start_time = time.time()

##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('/mnt/solar08b/c0287shimizu/R2D2_ss/run/d002/data/')
d.read_qq(10,'ro')    
d.read_qq(10,'op')    
d.read_qq_tau(100)  
d.read_qq(10, 'te')  
# d.read_qq_rt(10, 'in')

ro0 = d.p['ro0']      # background field density
ro1 = d.qq['ro']      # density disturbances
op = d.qq['op']       # opacity
x = d.p['x']          # height of index i
y = d.p['y']          # lenght of index j
z = d.p['z']          # lenght of index k
t1 = d.qq['te']       # temperature disturbances
t0 = d.p['te0']       # background of tenperature


#####  information #####
# MHD grid number
nx = 256
ny = 256
nz = 256
# RAD grid number
lx = nx - 1
ly = ny + 1
lz = nz + 1

# grid span
dxx = 1
dyy = 1
dzz = 1
dxxi = 1 / dxx
dyyi = 1 / dyy
dzzi = 1 / dzz

# ray direction
step = [-1, 1]                  # step direction
shift = [-1, 0]                 # shift of index
shift_s = [0, -1]               # not shift 
ista = [lx-1, 1]                # start point of i loop
iend = [1, lx-1]                # finish point of i loop
jsta = [ly-1, 1]                # start point of j loop
jend = [1, ly-1]                # finish point of j loop
ksta = [lz-1, 1]                # start point of k loop
kend = [1, lz-1]                # finish point of k loop

# short charcteristic
mrad = 3                        # number of rays in 1 octant
lr = np.zeros(mrad)             # lenght of ray in short characteristic
dc = np.zeros((3, 3))           # delta xc, yc (myu1(dxc, dyc, dzc), myu2(dxc, dyc, dzc), myu3(dxc, dyc, dzc))
dd = np.zeros((3, 3))           # delta xd, yd (myu1(dxd, dyd, dzd), myu2(dxd, dyd, dzd), myu3(dxd, dyd, dzd))
mux = np.zeros(mrad)            # x distance
muy = np.zeros(mrad)            # y distance
muz = np.zeros(mrad)            # z distance
omr = np.zeros(mrad)            # weight of ray

# ray distance
mux[0] = np.sqrt(7 / 9)
muy[0] = 1 / 3
muz[0] = 1 / 3

mux[1] = 1 / 3
muy[1] = np.sqrt(7 / 9)
muz[1] = 1 / 3

mux[2] = 1 / 3
muy[2] = 1 / 3
muz[2] = np.sqrt(7 / 9)

# ray weight
omr[:] = 1 / 24

# ray lenght
for m in range(mrad):
    lr[m] = min(dxx / mux[m], dyy / muy[m], dzz / muz[m])

# value to determine the index of each ray
dxc = np.zeros(mrad)
dyc = np.zeros(mrad)
dzc = np.zeros(mrad)
dxd = np.zeros(mrad)
dyd = np.zeros(mrad)
dzd = np.zeros(mrad)

di = np.zeros((mrad, 2, 4), dtype = int)
dj = np.zeros((mrad, 2, 4), dtype = int)
dk = np.zeros((mrad, 2, 4), dtype = int)

dc1 = np.zeros(mrad)
dd1 = np.zeros(mrad)
dc2 = np.zeros(mrad)
dd2 = np.zeros(mrad)
di1 = np.zeros(mrad)
di2 = np.zeros(mrad)

for m in range(mrad):
    dxc[m] = mux[m] * lr[m]
    dyc[m] = muy[m] * lr[m]
    dzc[m] = muz[m] * lr[m]

    dxd[m] = dxx - dxc[m]
    dyd[m] = dyy - dyc[m]
    dzd[m] = dzz - dzc[m]

    # y-z plane
    if dxx / mux[m] <= dyy / muy[m] and dxx / mux[m] <= dzz / muz[m]:
        dc1[m] = dyc[m]
        dd1[m] = dyd[m]
        dc2[m] = dzc[m]
        dd2[m] = dzd[m]
        di1[m] = 1.0 / dyy
        di2[m] = 1.0 / dzz
        di[m, 0, :] = 0
        di[m, 1, :] = -1
        dj[m, 0, :] = [ 0, -1,  0, -1]
        dk[m, 0, :] = [ 0,  0, -1, -1]
        dj[m, 1, :] = [-1,  0, -1,  0]
        dk[m, 1, :] = [-1, -1,  0,  0]

    # x-z plane
    elif dyy / muy[m] <= dxx / mux[m] and dyy / muy[m] <= dzz / muz[m]:
        dc1[m] = dxc[m]
        dd1[m] = dxd[m]
        dc2[m] = dzc[m]
        dd2[m] = dzd[m]
        di1[m] = 1.0 / dxx
        di2[m] = 1.0 / dzz
        dj[m, 0, :] = 0
        dj[m, 1, :] = -1
        di[m, 0, :] = [0, -1, 0, -1]
        dk[m, 0, :] = [0, 0, -1, -1]
        di[m, 1, :] = [-1, 0, -1, 0]
        dk[m, 1, :] = [-1, -1, 0, 0]

    # x-y plane
    elif dzz / muz[m] <= dxx / mux[m] and dzz / muz[m] <= dyy / muy[m]:
        dc1[m] = dxc[m]
        dd1[m] = dxd[m]
        dc2[m] = dyc[m]
        dd2[m] = dyd[m]
        di1[m] = 1.0 / dxx
        di2[m] = 1.0 / dyy
        dk[m, 0, :] = 0
        dk[m, 1, :] = -1
        di[m, 0, :] = [0, -1, 0, -1]
        dj[m, 0, :] = [0, 0, -1, -1]
        di[m, 1, :] = [-1, 0, -1, 0]
        dj[m, 1, :] = [-1, -1, 0, 0]

# list
alfa = np.zeros((lx, ly, lz))                                 # absorption_coefficient in RAD grid
B = np.zeros((lx, ly, lz))                                    # planck function in RAD grid
delta_tu = np.zeros((mrad, 2, 2, 2, lx, ly, lz))              # delta opacity
delta_I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))               # delta intensity
I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))                     # intensity
I_up = np.zeros((mrad, 2, 2, 2, lx, ly, lz))                  # intensity of upstream point
err_j = np.zeros((lx, lz))                                    # error of between jsta and jend
err_k = np.zeros((lx, ly))                                    # error of between ksta and kend

delta_tu_1ray = np.zeros((lx, ly, lz))                        # delta opacity for 1 ray
tu_1ray = np.zeros((lx, ly, lz))                              # opacity for 1 ray
x_tu1 = []                                                    # xindex of tu = 1
tu0 = 0.1                                
Fx = np.zeros((lx, ly, lz))                                   # flux in x direction
Fy = np.zeros((lx, ly, lz))                                   # flux in y direction
Fz = np.zeros((lx, ly, lz))                                   # flux in z direction
J = np.zeros((lx, ly, lz))                                    # mean intensity
Qrad = np.zeros((lx, ly, lz))                                 # radiative heating rate



##### define #####
# absorption cefficient for RAD grid
def absorption_coefficient_RAD(alfa_MHD, i, j, k):
    return exp((log(alfa_MHD[i,   j,   k]) + log(alfa_MHD[i,   j,   k+1]) \
              + log(alfa_MHD[i,   j+1, k]) + log(alfa_MHD[i,   j+1, k+1]) \
              + log(alfa_MHD[i+1, j,   k]) + log(alfa_MHD[i+1, j,   k+1]) \
              + log(alfa_MHD[i+1, j+1, k]) + log(alfa_MHD[i+1, j+1, k+1])) \
            / 8)

# planck function for RAD grid
def planck_function_RAD(B_MHD, i, j, k):
    return exp((log(B_MHD[i,   j,   k]) + log(B_MHD[i,   j,   k+1]) \
              + log(B_MHD[i,   j+1, k]) + log(B_MHD[i,   j+1, k+1]) \
              + log(B_MHD[i+1, j,   k]) + log(B_MHD[i+1, j,   k+1]) \
              + log(B_MHD[i+1, j+1, k]) + log(B_MHD[i+1, j+1, k+1])) \
            / 8)

# delta tu
def delta_opacity(alfa_up, alfa_down, log_alfa_up, log_alfa_down, lray):
    return ((alfa_down - alfa_up) * lray) / (log_alfa_down - log_alfa_up)

# delta intensity
def delta_intensity(B_up, B_down, log_B_up, log_B_down, delta_tu, m, idir, jdir, kdir, i, j, k):
    return ((B_down - B_up * exp(- delta_tu[m, idir, jdir, kdir, i, j, k])) * delta_tu[m, idir, jdir, kdir, i, j, k]) / (log_B_down - log_B_up + delta_tu[m, idir, jdir, kdir, i, j, k])

# Intensity at upstream
def intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2):
    return ((\
            I[m, idir, jdir, kdir, i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]] * dc1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]] * dd1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]] * dc1[m] * dd2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]] * dd1[m] * dd2[m] \
            )\
            * di1[m] * di2[m])



##### absorption coefficient in MHD grid #####
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
ro = Ro0 + ro1                                                        # density

alfa_data = ro * op                                                   # alfa_MHD_data is no Periodic Boundary Conditions 
alfa_MHD_ = np.insert(alfa_data, 0, alfa_data[:, ny-1, :], axis = 1)
alfa_MHD_ = np.insert(alfa_MHD_, ny+1, alfa_data[:, 0, :], axis = 1)
alfa_MHD = np.insert(alfa_MHD_, 0, alfa_MHD_[:, :, nz-1], axis = 2)   # alfa in MHD grid (Periodic Boundary Conditions)
alfa_MHD = np.insert(alfa_MHD, nz+1, alfa_MHD_[:, :, 0], axis = 2)

##### planck function in MHD grid#####
T0 , tmp, tmp = np.meshgrid(t0, d.p['y'], d.p['z'], indexing = 'ij')
t = T0 + t1                                                           # tempurature
sig = 5.67 * 10**(-5)                                                 # Stefan-Boltzmann constant

B_data = (sig * t**4) / pi                                            # B_MHD_data is no Periodic Boundary Conditions 
B_MHD_ = np.insert(B_data, 0, B_data[:, ny-1, :], axis = 1)
B_MHD_ = np.insert(B_MHD_, ny+1, B_data[:, 0, :], axis = 1)
B_MHD = np.insert(B_MHD_, 0, B_MHD_[:, :, nz-1], axis = 2)            # B in MHD grid (Periodic Boundary Conditions)
B_MHD = np.insert(B_MHD, nz+1, B_MHD_[:, :, 0], axis = 2)




###### absorption coefficient in RAD grid & planck function in RAD grid #####
print('##### alfa & planck function in RAD grid #####')
for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            alfa[i, j, k] = absorption_coefficient_RAD(alfa_MHD, i, j, k)
            B[i, j, k] = planck_function_RAD(B_MHD, i, j, k)




##### Opacity and Intensity #####
print('##### opacity and intensity #####')
### initial condition ###
# y-z plane
I[:, 0, :, :, lx-1, :, :] = 1e-10
I[:, 1, :, :, 0, :, :] = B[0, :, :]
# x-z plane
I[:, :, 0, :, :, ly-1, :] = B[:, ly-1, :]
I[:, :, 1, :, :, 0, :] = B[:, 0, :]
# x-y plane
I[:, :, :, 0, :, :, lz-1] = B[:, :, lz-1]
I[:, :, :, 1, :, :, 0] = B[:, :, 0]

for m in range(0, 1):                     
    for idir in range(0, 1):
        for jdir in range(0,1):
            for kdir in range(0, 1):
                
                print('### m, idir, jdir, kdir:', m, idir, jdir, kdir, ' ###')

                for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
                    for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                        for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                
                            ###### downstream #####
                            log_alfa_down = log(alfa[i + shift[idir], j + shift[jdir], k + shift[kdir]])    
                            alfa_down     =     alfa[i + shift[idir], j + shift[jdir], k + shift[kdir]]     
 
                            log_B_down = log(B[i + shift[idir], j + shift[jdir], k + shift[kdir]])          
                            B_down     =     B[i + shift[idir], j + shift[jdir], k + shift[kdir]]           


                            ##### upstream #####
                            log_alfa_up = ((\
                                           log(alfa[i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]]) * dc1[m] * dc2[m] \
                                         + log(alfa[i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]]) * dd1[m] * dc2[m] \
                                         + log(alfa[i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]]) * dc1[m] * dd2[m] \
                                         + log(alfa[i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]]) * dd1[m] * dd2[m] \
                                           ) * di1[m] * di2[m])
                
                            log_B_up = ((\
                                        log(B[i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]]) * dc1[m] * dc2[m] \
                                      + log(B[i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]]) * dd1[m] * dc2[m] \
                                      + log(B[i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]]) * dc1[m] * dd2[m] \
                                      + log(B[i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]]) * dd1[m] * dd2[m] \
                                        ) * di1[m] * di2[m])
                
                            alfa_up = exp(log_alfa_up)
                            B_up    = exp(log_B_up)
                            
                            # ray lenght lray[m] in short characteristic
                            if m == 0:
                                lray = lr[m] * (x[i] - x[i-1])
                            
                            elif m == 1:
                                lray = lr[m] * (y[i] - y[i-1])
                            
                            elif m == 2:
                                lray = lr[m] * (z[i] - z[i-1])

                            # delta opacity & delta intensity
                            delta_tu[m, idir, jdir, kdir, i, j, k] = delta_opacity(alfa_up, alfa_down, log_alfa_up, log_alfa_down, lray)
                            delta_I[m, idir, jdir, kdir, i, j, k] = delta_intensity(B_up, B_down, log_B_up, log_B_down, delta_tu, m, idir, jdir, kdir, i, j, k)


plot = delta_I[0, 0, 0, 0, 0, :, :]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

plt.savefig('test.png')
plt.close()