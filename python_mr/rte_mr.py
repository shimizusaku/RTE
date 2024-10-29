# This code is RTE for multi rays.
import os
import time
import R2D2
import numpy as np
# import matplotlib.pyplot as plt
from math import log, exp, pi, copysign

start_time = time.time()

##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('/mnt/solar07a/c0287shimizu/R2D2/run/d001/data/')
d.read_qq(10,'ro')    
d.read_qq(10,'op')    
d.read_qq(10, 'te')  
d.read_qq(10, 'se')  

ro0 = d.p['ro0']      # background field density
ro1 = d.qq['ro']      # density disturbances
op = d.qq['op']       # opacity
x = d.p['x']          # height of index i
y = d.p['y']          # lenght of index j
z = d.p['z']          # lenght of index k
t1 = d.qq['te']       # temperature disturbances
t0 = d.p['te0']       # background of tenperature
se0 = d.p['se0']      # background field entropy
se1 = d.qq['se']      # entropy disturbances


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
mux_1ray = 1 /sqrt(3)
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


# judge qrad
nx_jf_judge = 2
jf_judge = [1, 0]
tumin = 1e-2
tumax = 2e-1
log_tumin = log(tumin)
log_tumax = log(tumax)
delta_log_tu = (log_tumax - log_tumin) / (nx_jf_judge - 1)
jk_coef = delta_log_tu / (log_tumin - log_tumax)


# list
alpha = np.zeros((lx, ly, lz))                                 # absorption_coefficient in RAD grid
beta = np.zeros((lx, ly, lz))                                    # planck function in RAD grid
delta_tu = np.zeros((mrad, 2, 2, 2, lx, ly, lz))              # delta opacity
exp_delta_tu = np.zeros((mrad, 2, 2, 2, lx, ly, lz))          # exp(delta_tu)
delta_I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))               # delta intensity
I = np.zeros((mrad, 2, 2, 2, lx, ly, lz))                     # intensity
I_up = np.zeros((mrad, 2, 2, 2, lx, ly, lz))                  # intensity of upstream point
err_j = np.zeros((lx, lz))                                    # error of between jsta and jend
err_k = np.zeros((lx, ly))                                    # error of between ksta and kend

delta_tu_1ray = np.zeros((lx, ly, lz))                        # delta opacity for 1 ray
tu_1ray = np.zeros((nx, ly, lz))                              # opacity for 1 ray
x_tu1 = np.zeros((ly, lz))                                                    # xindex of tu = 1i
Qrad_tu1 = np.zeros((ly, lz))
tu0 = 0.1                                
Fx = np.zeros((lx, ly, lz))                                   # flux in x direction
Fy = np.zeros((lx, ly, lz))                                   # flux in y direction
Fz = np.zeros((lx, ly, lz))                                   # flux in z direction
J = np.zeros((lx, ly, lz))                                    # mean intensity
Qrad = np.zeros((lx, ly, lz))                                 # radiative heating rate



##### define #####
# absorption cefficient for RAD grid
def absorption_coefficient_RAD(alpha_MHD, i, j, k):
    return exp((log(alpha_MHD[i,   j,   k]) + log(alpha_MHD[i,   j,   k+1]) \
              + log(alpha_MHD[i,   j+1, k]) + log(alpha_MHD[i,   j+1, k+1]) \
              + log(alpha_MHD[i+1, j,   k]) + log(alpha_MHD[i+1, j,   k+1]) \
              + log(alpha_MHD[i+1, j+1, k]) + log(alpha_MHD[i+1, j+1, k+1])) \
            / 8)

# planck function for RAD grid
def planck_function_RAD(beta_mhd, i, j, k):
    return exp((log(beta_mhd[i,   j,   k]) + log(beta_mhd[i,   j,   k+1]) \
              + log(beta_mhd[i,   j+1, k]) + log(beta_mhd[i,   j+1, k+1]) \
              + log(beta_mhd[i+1, j,   k]) + log(beta_mhd[i+1, j,   k+1]) \
              + log(beta_mhd[i+1, j+1, k]) + log(beta_mhd[i+1, j+1, k+1])) \
            / 8)

# Intensity at upstream
def intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2):
    return ((\
            I[m, idir, jdir, kdir, i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]] * dc1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]] * dd1[m] * dc2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]] * dc1[m] * dd2[m] \
          + I[m, idir, jdir, kdir, i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]] * dd1[m] * dd2[m] \
            )\
            * di1[m] * di2[m])


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
def rte_delta_I(beta_down, beta_up, log_beta_down, log_beta_up, \
                 delta_tu, m, idir, jdir, kdir, i, j, k):
    # parameter to check
    epc = 1e-4

    exp_delta_tu = rte_exp(delta_tu, m, idir, jdir, kdir, i, j, k)

    beta_upt = beta_up * exp_delta_tu

    beta_du = beta_upt / beta_down - 1
    beta_dum = copysign(1, beta_du) * max(abs(beta_du), epc) 

    beta_du15 = 1 - 0.5 * beta_du
    # for zero divide
    beta_du15m = copysign(1, beta_du15) * max(abs(beta_du15), epc) 

    tmp = 0.5 + copysign(0.5, abs(beta_du) - epc)
    log_beta_du_delta_tu0 = log_beta_down - log_beta_up + delta_tu[m, idir, jdir, kdir, i, j, k]
    log_beta_du_delta_tu1 = 1                                       # value is not important

    log_beta_du_delta_tu = \
        log_beta_du_delta_tu0 * tmp + log_beta_du_delta_tu1 * (1 - tmp)

    delta_in0 = (beta_down - beta_upt) / log_beta_du_delta_tu
    delta_in1 = beta_down / beta_du15m

    delta_in = (delta_in0 * tmp + delta_in1 * (1 - tmp)) * delta_tu[m, idir, jdir, kdir, i, j, k]
    
    return delta_in

# exp(x)
def rte_exp(delta_tu, m, idir, jdir, kdir, i, j, k):
    xmin = -700
    tmp = 0.5 + copysign(0.5, - delta_tu[m, idir, jdir, kdir, i, j, k] - xmin)

    return exp(max(- delta_tu[m, idir, jdir, kdir, i, j, k], xmin)) * tmp



##### absorption coefficient in MHD grid #####
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
ro = Ro0 + ro1                                                        # density

alpha_data = ro * op                                                   # alpha_MHD_data is no Periodic Boundary Conditions 
alpha_MHD_ = np.insert(alpha_data, 0, alpha_data[:, ny-1, :], axis = 1)
alpha_MHD_ = np.insert(alpha_MHD_, ny+1, alpha_data[:, 0, :], axis = 1)
alpha_MHD = np.insert(alpha_MHD_, 0, alpha_MHD_[:, :, nz-1], axis = 2)   # alpha in MHD grid (Periodic Boundary Conditions)
alpha_MHD = np.insert(alpha_MHD, nz+1, alpha_MHD_[:, :, 0], axis = 2)


##### planck function in MHD grid#####
T0 , tmp, tmp = np.meshgrid(t0, d.p['y'], d.p['z'], indexing = 'ij')
t = T0 + t1                                                           # tempurature
sig = 5.67 * 10**(-5)                                                 # Stefan-Boltzmann constant

beta_data = (sig * t**4) / pi                                            # beta_mhd_data is no Periodic Boundary Conditions 
beta_mhd_ = np.insert(beta_data, 0, beta_data[:, ny-1, :], axis = 1)
beta_mhd_ = np.insert(beta_mhd_, ny+1, beta_data[:, 0, :], axis = 1)
beta_mhd = np.insert(beta_mhd_, 0, beta_mhd_[:, :, nz-1], axis = 2)            # B in MHD grid (Periodic Boundary Conditions)
beta_mhd = np.insert(beta_mhd, nz+1, beta_mhd_[:, :, 0], axis = 2)


##### entropy #####
Se0, tmp, tmp = np.meshgrid(se0, d.p['y'], d.p['z'], indexing='ij')
se = Se0 + se1                                                         # entropy

se_mhd_ = np.insert(se, 0, se[:, ny-1, :], axis = 1)
se_mhd_ = np.insert(se_mhd_, ny+1, se[:, 0, :], axis = 1)
se_mhd = np.insert(se_mhd_, 0, se_mhd_[:, :, nz-1], axis = 2)            
se_mhd = np.insert(se_mhd, nz+1, se_mhd_[:, :, 0], axis = 2)


###### absorption coefficient in RAD grid & planck function in RAD grid #####
print('##### alpha & planck function in RAD grid #####')
for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            alpha[i, j, k] = absorption_coefficient_RAD(alpha_MHD, i, j, k)
            beta[i, j, k] = planck_function_RAD(beta_mhd, i, j, k)




##### Opacity and Intensity #####
print('##### opacity and intensity #####')
### initial condition ###
# y-z plane
I[:, 0, :, :, lx-1, :, :] = 1e-10
I[:, 1, :, :, 0, :, :] = beta[0, :, :]
# x-z plane
I[:, :, 0, :, :, ly-1, :] = beta[:, ly-1, :]
I[:, :, 1, :, :, 0, :] = beta[:, 0, :]
# x-y plane
I[:, :, :, 0, :, :, lz-1] = beta[:, :, lz-1]
I[:, :, :, 1, :, :, 0] = beta[:, :, 0]

for m in range(0, 3):                     
    for idir in range(0, 2):
        for jdir in range(0, 2):
            for kdir in range(0, 2):
                
                print('### m, idir, jdir, kdir:', m, idir, jdir, kdir, ' ###')

                for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
                    for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                        for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                
                            ###### downstream #####
                            log_alpha_down = log(alpha[i + shift[idir], j + shift[jdir], k + shift[kdir]])    
                            alpha_down     =     alpha[i + shift[idir], j + shift[jdir], k + shift[kdir]]     
 
                            log_beta_down = log(beta[i + shift[idir], j + shift[jdir], k + shift[kdir]])          
                            beta_down     =     beta[i + shift[idir], j + shift[jdir], k + shift[kdir]]           


                            ##### upstream #####
                            log_alpha_up = ((\
                                           log(alpha[i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]]) * dc1[m] * dc2[m] \
                                         + log(alpha[i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]]) * dd1[m] * dc2[m] \
                                         + log(alpha[i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]]) * dc1[m] * dd2[m] \
                                         + log(alpha[i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]]) * dd1[m] * dd2[m] \
                                           ) * di1[m] * di2[m])
                
                            log_beta_up = ((\
                                        log(beta[i + di[m, idir, 0], j + dj[m, jdir, 0], k + dk[m, kdir, 0]]) * dc1[m] * dc2[m] \
                                      + log(beta[i + di[m, idir, 1], j + dj[m, jdir, 1], k + dk[m, kdir, 1]]) * dd1[m] * dc2[m] \
                                      + log(beta[i + di[m, idir, 2], j + dj[m, jdir, 2], k + dk[m, kdir, 2]]) * dc1[m] * dd2[m] \
                                      + log(beta[i + di[m, idir, 3], j + dj[m, jdir, 3], k + dk[m, kdir, 3]]) * dd1[m] * dd2[m] \
                                        ) * di1[m] * di2[m])
                
                            alpha_up = exp(log_alpha_up)
                            beta_up    = exp(log_beta_up)
                            
                            # ray lenght lray[m] in short characteristic
                            if m == 0:
                                lray = lr[m] * (x[i] - x[i-1])
                            
                            elif m == 1:
                                lray = lr[m] * (y[i] - y[i-1])
                            
                            elif m == 2:
                                lray = lr[m] * (z[i] - z[i-1])

                            # delta opacity & delta intensity
                            # delta_tu[m, idir, jdir, kdir, i, j, k] = delta_opacity(alpha_up, alpha_down, log_alpha_up, log_alpha_down, lray)
                            delta_tu[m, idir, jdir, kdir, i, j, k] = rte_delta_tu(alpha_down, alpha_up, lray)
                            exp_delta_tu[m, idir, jdir, kdir, i, j, k] = rte_exp(delta_tu, m, idir, jdir, kdir, i, j, k)
                            # delta_I[m, idir, jdir, kdir, i, j, k] = delta_intensity(beta_up, beta_down, log_beta_up, log_beta_down, delta_tu, m, idir, jdir, kdir, i, j, k, exp_delta_tu)
                            delta_I[m, idir, jdir, kdir, i, j, k] = rte_delta_I(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu, m, idir, jdir, kdir, i, j, k)


                count = 1   # The count refers to the number of iterations in the calculation.

                if m == 0:
                    converged = False
                    while not converged:                    
                        for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
                            for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp_delta_tu[m, idir, jdir, kdir, i, j, k] + delta_I[m, idir, jdir, kdir, i, j, k]
        
                        # error between sta and end
                        for i in range(1 + shift[idir], lx + shift[idir]):
                            for j in range(0, ly):
                                for k in range(0, lz):
                                    Irecv_j = I[m, idir, jdir, kdir, i, jsta[jdir] + shift_s[jdir], k] 
                                    Isend_j = I[m, idir, jdir, kdir, i, jend[jdir] + shift[jdir], k]
                                    err_j[i, k] = abs(Isend_j - Irecv_j) / Isend_j        

                                    Irecv_k = I[m, idir, jdir, kdir, i, j, ksta[kdir] + shift_s[kdir]]
                                    Isend_k = I[m, idir, jdir, kdir, i, j, kend[kdir] + shift[kdir]  ] 
                                    err_k[i, j] = abs(Isend_k - Irecv_k) / Isend_k

                        max_err_j = np.amax(err_j[:, :])
                        max_err_k = np.amax(err_k[:, :])
                        max_err = max(max_err_j, max_err_k)    

                        if max_err < 10**(-3):                            
                            # Terminate the iterative calculation.
                            print('count:', count)
                            converged = True

                        else:
                            # update of sta value by end value
                            I[m, idir, jdir, kdir, :, jsta[jdir] + shift_s[jdir], :] = I[m, idir, jdir, kdir, :, jend[jdir] + shift[jdir], :]
                            I[m, idir, jdir, kdir, :, :, ksta[kdir] + shift_s[kdir]] = I[m, idir, jdir, kdir, :, :, kend[kdir] + shift[kdir]] 
                            
                            count = count + 1 


                elif m == 1:
                    converged = False
                    while not converged:
                        for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):                                            
                            for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
                                for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):
                                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp_delta_tu[m, idir, jdir, kdir, i, j, k] + delta_I[m, idir, jdir, kdir, i, j, k]
        
                        # error between sta and end
                        for i in range(1 + shift[idir], lx + shift[idir]):
                            for j in range(0, ly):
                                for k in range(0, lz):
                                    Irecv_j = I[m, idir, jdir, kdir, i, jsta[jdir] + shift_s[jdir], k] 
                                    Isend_j = I[m, idir, jdir, kdir, i, jend[jdir] + shift[jdir],   k] 
                                    err_j[i, k] = (Isend_j - Irecv_j) / Isend_j        

                                    Irecv_k = I[m, idir, jdir, kdir, i, j, ksta[kdir] + shift_s[kdir]]
                                    Isend_k = I[m, idir, jdir, kdir, i, j, kend[kdir] + shift[kdir]  ] 
                                    err_k[i, j] = (Isend_k - Irecv_k) / Isend_k

                        max_err_j = np.amax(abs(err_j[:, :]))
                        max_err_k = np.amax(abs(err_k[:, :]))
                        max_err = max(max_err_j, max_err_k)

                        if max_err < 10**(-3):
                            # Terminate the iterative calculation.
                            print('count:', count)
                            converged = True

                        else:
                            # update of sta value by end value
                            I[m, idir, jdir, kdir, :, jsta[jdir] + shift_s[jdir], :] = I[m, idir, jdir, kdir, :, jend[jdir] + shift[jdir], :]
                            I[m, idir, jdir, kdir, :, :, ksta[kdir] + shift_s[kdir]] = I[m, idir, jdir, kdir, :, :, kend[kdir] + shift[kdir]] 

                            count = count + 1 


                elif m == 2:
                    converged = False
                    while not converged:
                        for k in range(ksta[kdir], kend[kdir] + step[kdir], step[kdir]):                                           
                            for i in range(ista[idir], iend[idir] + step[idir], step[idir]):
                                for j in range(jsta[jdir], jend[jdir] + step[jdir], step[jdir]):
                                    I_up[m, idir, jdir, kdir, i, j, k] = intensity_up(I, m, idir, jdir, kdir, di, dj, dk, i, j, k, dc1, dc2, dd1, dd2, di1, di2)
                                    I[m, idir, jdir, kdir, i + shift[idir], j + shift[jdir], k + shift[kdir]] = I_up[m, idir, jdir, kdir, i, j, k] * exp_delta_tu[m, idir, jdir, kdir, i, j, k] + delta_I[m, idir, jdir, kdir, i, j, k]
        
                        # error between sta and end
                        for i in range(1 + shift[idir], lx + shift[idir]):
                            for j in range(0, ly):
                                for k in range(0, lz):
                                    Irecv_j = I[m, idir, jdir, kdir, i, jsta[jdir] + shift_s[jdir], k] 
                                    Isend_j = I[m, idir, jdir, kdir, i, jend[jdir] + shift[jdir],   k] 
                                    err_j[i, k] = (Isend_j - Irecv_j) / Isend_j      

                                    Irecv_k = I[m, idir, jdir, kdir, i, j, ksta[kdir] + shift_s[kdir]]
                                    Isend_k = I[m, idir, jdir, kdir, i, j, kend[kdir] + shift[kdir]  ] 
                                    err_k[i, j] = (Isend_k - Irecv_k) / Isend_k

                        max_err_j = np.amax(abs(err_j[:, :]))
                        max_err_k = np.amax(abs(err_k[:, :]))
                        max_err = max(max_err_j, max_err_k)

                        if max_err < 10**(-3):
                            # Terminate the iterative calculation.
                            print('count:', count)
                            converged = True

                        else:
                            # update of sta value by end value
                            I[m, idir, jdir, kdir, :, jsta[jdir] + shift_s[jdir], :] = I[m, idir, jdir, kdir, :, jend[jdir] + shift[jdir], :]
                            I[m, idir, jdir, kdir, :, :, ksta[kdir] + shift_s[kdir]] = I[m, idir, jdir, kdir, :, :, kend[kdir] + shift[kdir]] 

                            count = count + 1 



##### radiative heating rate Q #####
print('##### radiative heating rate #####')
for i in range(lx-1, -1, -1):
    for j in range(1, ly):
        for k in range(1, lz):
            # downstream
            alpha_down = alpha_MHD[i, j, k]

            # upstream
            alpha_up = alpha_MHD[i+1, j, k]

            # ray lenght lray(m) in short characteristics
            lray = (x[i+1] - x[i]) / mux_1ray

            # delta optical depth 1ray & optical depth 1ray
            delta_tu_1ray[i, j, k] = rte_delta_tu(alpha_down, alpha_up, lray)
            # delta_tu_1ray[i, j, k] = ((alpha_MHD[i, j, k] - alpha_MHD[i-1, j, k]) * (x[i] - x[i-1])) / (log(alpha_MHD[i, j, k]) - log(alpha_MHD[i-1, j, k]))

            # optical depth 1ray
            tu_1ray[i, j, k] = tu_1ray[i+1, j, k] + delta_tu_1ray[i, j, k]

            if (1 - tu_1ray[i, j, k]) * (1 - tu_1ray[i+1, j, k]) <= 0:
                x_tu1[j, k] = i




for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            for m in range(0, 3):
                for idir in range(0, 2):
                    for jdir in range(0, 2):
                        for kdir in range(0, 2):
                            Fx[i, j, k] = Fx[i, j, k] + (4 * pi * omr[m] * mux[m] * step[idir] * I[m, idir, jdir, kdir, i, j, k])
                            Fy[i, j, k] = Fy[i, j, k] + (4 * pi * omr[m] * muy[m] * step[jdir] * I[m, idir, jdir, kdir, i, j, k])
                            Fz[i, j, k] = Fz[i, j, k] + (4 * pi * omr[m] * muz[m] * step[kdir] * I[m, idir, jdir, kdir, i, j, k])

                            J[i, j, k] = J[i, j, k] + (omr[m] * I[m, idir, jdir, kdir, i, j, k])

for i in range(1, lx):
    for j in range(1, ly):
        for k in range(1, lz):
            # Q_F 
            # up is value in upstream point, down is value in downstream point)
            Fx_up   = (Fx[i,   j, k] + Fx[i,   j-1, k] + Fx[i,   j, k-1] + Fx[i,   j-1, k-1]) / 4
            Fx_down = (Fx[i-1, j, k] + Fx[i-1, j-1, k] + Fx[i-1, j, k-1] + Fx[i-1, j-1, k-1]) / 4

            Fy_up   = (Fy[i, j,   k] + Fy[i-1, j,   k] + Fy[i, j,   k-1] + Fy[i-1, j,   k-1]) / 4
            Fy_down = (Fy[i, j-1, k] + Fy[i-1, j-1, k] + Fy[i, j-1, k-1] + Fy[i-1, j-1, k-1]) / 4

            Fz_up   = (Fz[i, j, k  ] + Fz[i-1, j, k  ] + Fz[i, j-1, k  ] + Fz[i-1, j-1, k  ]) / 4
            Fz_down = (Fz[i, j, k-1] + Fz[i-1, j, k-1] + Fz[i, j-1, k-1] + Fz[i-1, j-1, k-1]) / 4

            Q_F = - (Fx_up - Fx_down) / (x[i] - x[i-1]) \
                  - (Fy_up - Fy_down) / (y[i] - y[i-1]) \
                  - (Fz_up - Fz_down) / (z[i] - z[i-1])
            
            # Q_J
            J_MHD = (J[i,   j, k] + J[i,   j, k-1] + J[i,   j-1, k] + J[i,   j-1, k-1] \
                   + J[i-1, j, k] + J[i-1, j, k-1] + J[i-1, j-1, k] + J[i-1, j-1, k-1]) / 8
            
            Q_J = 4 * pi * alpha_data[i, j-1, k-1] * (J_MHD - beta_data[i, j-1, k-1])

            # Qrad
            # tu_norm = exp(- (tu_1ray[i, j, k] / tu0))
            # Qrad[i, j, k] = (Q_J * tu_norm + Q_F * (1 - tu_norm)) / (ro[i, j-1, k-1] * t[i, j-1, k-1])
            
            tmp = min(max((log(tu_1ray[i, j, k]) - log_tumin) / delta_log_tu + 1, 1), nx_jf_judge)
            tmp = (tmp - int(tmp)) * jk_coef + jf_judge[int(tmp) - 1]
            
            Qrad[i, j, k] = ((Q_J * tmp + Q_F * (1 - tmp)) / (ro[i, j-1, k-1] * t[i, j-1, k-1])) \
                          * (0.5 + copysign(0.5, tu_1ray[i, j, k] - 1e-4)) \
                          + (se_mhd[i, j-1, k-1] / t[i, j-1, k-1]) * max(0, (2000 - t[i, j-1, k-1]/1))


# Qrad in the plane of tu=1
for j in range(1, ly):
    for k in range(1, lz):
        Qrad_tu1[j, k] = Qrad[int(x_tu1[j, k]), j, k]


##### save data #####
print('### save data ###')
np.savez('/mnt/solar07a/c0287shimizu/work/data/d001/MultiRay/python.npz', qrad=Qrad, intensity=I, tu_1ray=tu_1ray, xtu1=x_tu1, qrad_tu1_mr=Qrad_tu1)


##### plot #####
print('### plot ###')

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
                plot = I[m, idir, jdir, kdir, x_point[idir], :, :]

                plt.rcParams['xtick.direction'] = 'in'
                plt.rcParams['ytick.direction'] = 'in'

                plt.figure()
                plt.imshow(plot)
                plt.title(mti[m] + iti[idir] + jti[jdir] + kti[kdir])

                ticks = np.linspace(plot.min(), plot.max(), 10)
                cbar = plt.colorbar()
                cbar.set_ticks(ticks)

                save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_py/MultiRay'
                file_name = mti[m] + iti[idir] + jti[jdir] + kti[kdir] + '.png'
                file_path = os.path.join(save_dir, file_name)

                plt.savefig(file_path)
                plt.close()


# Qrad
plot = Qrad_tu1[1:, 1:]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(plot)
plt.title('Qrad tu1')

ticks = np.linspace(plot.min(), plot.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/mnt/solar07a/c0287shimizu/work/fig/d001_py/MultiRay'
file_name = 'Qrad_tu1.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()


end_time = time.time()
elapsed_time = end_time - start_time
print('Calculation time:', elapsed_time)