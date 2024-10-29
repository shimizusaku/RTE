# multi rays RTE
# FOA
# change code

import os
import time
import R2D2
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, pi, sqrt

start_time = time.time()

##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('../../run/d005/data/')
d.read_qq(10,'ro')   # data of density in 10 step
d.read_qq(10,'op')   # data of opasity in 10 step
d.read_qq_tau(600)   # ifac = 60

ro0 = d.p['ro0']      # background field density
ro1 = d.qq['ro']      # density disturbances
op = d.qq['op']       # opacity
# x = d.p['x']          # height of index i

d.read_qq(10, 'te')   # data of temperature in 10 step
t1 = d.qq['te']       # temperature disturbances
t0 = d.p['te0']       # background of tenperature



### grid information ###
# delta 
D = [1, 1, 1]            # cell width (dx, dy, dz)
dl = np.zeros(3)         # delta lenght (l_myu1, l_myu2, l_myu3)
dc = np.zeros((3, 3))    # delta xc, yc (myu1(dxc, dyc, dzc), myu2(dxc, dyc, dzc), myu3(dxc, dyc, dzc))
dd = np.zeros((3, 3))    # delta xd, yd (myu1(dxd, dyd, dzd), myu2(dxd, dyd, dzd), myu3(dxd, dyd, dzd))


# MHD grid number
nx = 64
ny = 64
nz = 64
# RAD grid number
lx = nx - 1
ly = ny + 1
lz = nz + 1

# ray information
mrad = 3                 # number of rays in 1 octant
step = [-1, 1]           # step direction
shift = [-1, 0]          # shift of index
ista = [lx-1, 0]         # start point of i roop
iend = [0, lx-1]         # finish point of i loop
jsta = [ly-1, 0]         # start point of j loop
jend = [0, ly-1]         # finish point of j loop
ksta = [lz-1, 0]         # start point of k loop
kend = [0, lz-1]         # finish point of k loop




### grid span for ray ###
mux = np.zeros(mrad)
muy = np.zeros(mrad)
muz = np.zeros(mrad)
omr = np.zeros(mrad)
lr = np.zeros(mrad)

# grid span
dxx = 1
dyy = 1
dzz = 1
dxxi = 1 / dxx
dyyi = 1 / dyy
dzzi = 1 / dzz

# ray direction
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

# span for each direction
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
        dj[m, 0, :] = [0, -1, 0, -1]
        dk[m, 0, :] = [0, 0, -1, -1]
        dj[m, 1, :] = [-1, 0, -1, 0]
        dk[m, 1, :] = [-1, -1, 0, 0]

    # x-z plane
    elif dyy / muy[m] <= dxx / mux[m] and dyy / muy[m] <= dzz / muz[m]:
        dc1[m] = dxc[m]
        dd1[m] = dxd[m]
        dc2[m] = dyc[m]
        dd2[m] = dyd[m]
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
def delta_opacity(alfa_up, lray):
    return alfa_up * lray

# delta intensity
def delta_intensity(B_down, delta_tu, m, i, j, k, idir, kdir, jdir):
    return B_down * (1 - exp(- delta_tu[m, idir, jdir, kdir, i, j, k]))

# Intensity at upstream
def intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2):
    return ((\
            I[m, idir, jdir, kdir, i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]] * dc1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]] * dd1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]] * dc1[m] * dd2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]] * dd1[m] * dd2[m] \
            )\
            * di1[m] * di2[m])




##### opacity #####
# density
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
ro = Ro0 + ro1   

# absorption_coefficient_MHD 
alfa_data = ro * op     # alfa_MHD_data is no Periodic Boundary Conditions 

alfa_MHD_ = np.insert(alfa_data, 0, alfa_data[:, ny-1, :], axis = 1)
alfa_MHD_ = np.insert(alfa_MHD_, ny+1, alfa_data[:, 0, :], axis = 1)

alfa_MHD = np.insert(alfa_MHD_, 0, alfa_MHD_[:, :, nz-1], axis = 2)   # alfa in MHD grid (Periodic Boundary Conditions)
alfa_MHD = np.insert(alfa_MHD, nz+1, alfa_MHD_[:, :, 0], axis = 2)

# temperature
T0 , tmp, tmp = np.meshgrid(t0, d.p['y'], d.p['z'], indexing = 'ij')
t = T0 + t1

# planck function
sig = 5.67 * 10**(-5)            # Stefan-Boltzmann constant [erg cm^-2 deg^-4 s^-1]
B_data = (sig * t**4) / pi       # B_MHD_data is no Periodic Boundary Conditions 

B_MHD_ = np.insert(B_data, 0, B_data[:, ny-1, :], axis = 1)
B_MHD_ = np.insert(B_MHD_, ny+1, B_data[:, 0, :], axis = 1)

B_MHD = np.insert(B_MHD_, 0, B_MHD_[:, :, nz-1], axis = 2)   # B in MHD grid (Periodic Boundary Conditions)
B_MHD = np.insert(B_MHD, nz+1, B_MHD_[:, :, 0], axis = 2)


# absorption coefficient RAD & planck function RAD
print('### alfa & planck function for RAD point ###')

alfa = np.zeros((lx, ly, lz))           # for absorption_coefficient RAD grid
B = np.zeros((lx, ly, lz))              # for planck function RAD grid

for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            alfa[i, j, k] = absorption_coefficient_RAD(alfa_MHD, i, j, k)
            B[i, j, k] = planck_function_RAD(B_MHD, i, j, k)


##### ray direction #####
delta_tu = np.zeros((mrad, 2, 2, 2, lx, ly, lz))
delta_I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))
I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))
I_up = np.zeros((mrad, 2, 2, 2, lx, ly, lz))
I_up_new = np.zeros((mrad, 2, 2, 2, lx, ly, lz))

log_alfa_up = np.zeros((lx, ly, lz))
C = []

### initial condition ###
lxb = 50

# y-z plane
I[:, 0, :, :, lx-1, :, :] = 0
I[:, 1, :, :, lxb, :, :] = B[lxb, :, :]

# x-z plane
I[1, :, :, :, :, ly-1, :] = B[:, ly-1, :]
I[1, :, :, :, :, 0, :] = B[:, 0, :]

# x-y plane
I[2, :, :, :, :, :, lz-1] = B[:, :, lz-1]
I[2, :, :, :, :, :, 0] = B[:, :, 0]


ista_I = [lx-1, lxb]              # start point of i loop for intensity
iend_I = [lxb, lx-1]              # finish point of i loop for intensity


for mray in range(0, 8 * mrad):                     # number of all rays 0, 8 * mrad
    print('ray number :', mray)
    m = mray // 8                                   # index for number of ray set(myu_m)
    idir = (mray - m * 8) // 4                      # index for x direction
    jdir = (mray - m * 8 - idir * 4) // 2           # index for y direction
    kdir = (mray - m * 8 - idir * 4 - jdir * 2)     # index for z direction

    # print(idir, jdir, kdir)

    count = 1
    for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
        for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
            for k in range(ksta[kdir], kend[kdir] + 1, step[kdir]):
                
                ###### downstream #####
                log_alfa_down = log(alfa[i + shift[idir], j + shift[jdir], k + shift[kdir]])    # log(alfa) in downstream
                alfa_down = alfa[i + shift[idir], j + shift[jdir], k + shift[kdir]]             # alfa in downstream
 
                log_B_down = log(B[i + shift[idir], j + shift[jdir], k + shift[kdir]])          # log(B) in downstream
                B_down = B[i + shift[idir], j + shift[jdir], k + shift[kdir]]                   # B in downstream

                # print(mray, i, j, k, alfa_down, B_down)

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

                # a[m, idir, jdir, kdir, i, j, k] = alfa_up
                B_up = exp(log_B_up)

                # print(mray, i, j, k, log_alfa_up, log_B_up)

                lray = lr[m] * 4800000

                delta_tu[m, idir, jdir, kdir, i, j, k] = delta_opacity(alfa_up, lray)
                # print('delta_tu:', delta_tu[m, idir, jdir, kdir, i, j, k])

                delta_I[m, idir, jdir, kdir, i, j, k] = delta_intensity(B_down, delta_tu, m, i, j, k, idir, kdir, jdir)
                # print(m, idir, jdir, kdir, i, j, k, delta_I[m, idir, jdir, kdir, i, j, k])

    # mrad = 1
    if m == 0:
        for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
            for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]

                    # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])
        
        print('### end ###')

   # mrad = 2
    elif m == 1:
        for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
            for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]
                
                    # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])

        # error calculation and updating I_up
        converged = False
        while not converged:
            for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                    for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                        I_up_new[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                    
            err = I_up[m, idir, jdir, kdir, :, jsta[jdir] + dj[m, jdir, 0], :] - I_up_new[m, idir, jdir, kdir, :, jsta[jdir] + dj[m, jdir, 0], :]
            max_err = np.amax(err)

            if max_err < 10**(-3):
                I_up[m, idir, jdir, kdir, :, :, :] = I_up_new[m, idir, jdir, kdir, :, :, :]

                for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                    for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                        for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                            I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]
                
                            # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])

                C.append(count)
                print('count:', count)
                print('### end ###')
                converged = True
                        
            else:
                count = count + 1
                
                I_up[1, idir, jdir, kdir, :, :, :] = I_up_new[1, idir, jdir, kdir, :, :, :]

                for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                    for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                        for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                            I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]
                            # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])
    
    # mrad = 3
    elif m == 2:
        for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
            for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]

                
                    # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])

        # error calculation and updating I_up
        converged = False
        while not converged:
            for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                    for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                        I_up_new[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                    
            err = I_up[m, idir, jdir, kdir, :, :, ksta[kdir]] - I_up_new[m, idir, jdir, kdir, :, :, ksta[kdir]]
            max_err = np.amax(err)

            if max_err < 10**(-3):
                I_up[m, idir, jdir, kdir, :, :, :] = I_up_new[m, idir, jdir, kdir, :, :, :]

                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                    for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                        for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                            I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]
                
                            # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])

                C.append(count)
                print('count:', count)
                print('### end ###')
                converged = True

            else:
                count = count + 1
                
                I_up[m, idir, jdir, kdir, :, :, :] = I_up_new[m, idir, jdir, kdir, :, :, :]

                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                    for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                        for i in range(ista_I[idir], iend_I[idir] + step[idir], step[idir]):
                            I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp(- delta_tu[m, idir, jdir, kdir, i, j, k]) + delta_I[m, idir, jdir, kdir, i, j, k]
                            # print(m, i + shift[idir], j + shift[jdir], k + shift[kdir], I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]])


##### radiative heatind rate Q #####
print('### radiative heating rate Q ###')




##### plot #####
print('### plot ###')
Intensity = I[2, 1, 1, 1, 62, :, :]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(Intensity)

ticks = np.linspace(Intensity.min(), Intensity.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/cidashome/sc/c0287shimizu/R2D2_ex/py/fig/MR/MR_test'
file_name = 'test.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()


end_time = time.time()

elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")