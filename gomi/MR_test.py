# multi rays RTE
# FOA
# test file (3 rays only)

import os
import R2D2
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, pi, sqrt



##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('../../run/d004/data/')
d.read_qq(10,'ro')   # data of density in 10 step
d.read_qq(10,'op')   # data of opasity in 10 step
d.read_qq_tau(600)   # ifac = 60

ro0 = d.p['ro0']      # background field density
ro1 = d.qq['ro']      # density disturbances
op = d.qq['op']       # opacity
x = d.p['x']          # height of index i

d.read_qq(10, 'te')   # data of temperature in 10 step
t1 = d.qq['te']       # temperature disturbances
t0 = d.p['te0']       # background of tenperature


### Carlson+1963 A4 integral ###
myu = [[sqrt(7/9), 1/3, 1/3], [1/3, sqrt(7/9), 1/3], [1/3, 1/3, sqrt(7/9)]]


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


# list
alfa = np.zeros((lx, ly, lz))           # for absorption_coefficient_RAD
alfa_up = np.zeros((3, lx-1, ly, lz))   # for absorption_coefficient_up (myu_num, i, j, k)
delta_tu = np.zeros((3, lx-1, ly, lz))  # for delta opacity
tu = np.zeros((3, lx, ly, lz))          # opacity
tu_up = np.zeros((3, lx-1, ly, lz))     # opacity for upwind point

B = np.zeros((lx, ly, lz))
B_up = np.zeros((3, lx-1, ly, lz))      # for planck_function_up 
delta_I = np.zeros((3, lx-1, ly, lz))   # for delta intensity
I = np.zeros((3, lx, ly, lz))           # intensity
I_up = np.zeros((3, lx-1, ly, lz))      # intensity for upwind point
I_up_new = np.zeros((3, lx-1, ly, lz))  # for check the error

# J_RAD = np.zeros((lx, ly, lz))          # mean intensity for RAD grid
J_cen = np.zeros((nx, ny, nz))          # mean intensity for cell center
Q_J = np.zeros((nx, ny, nz))            # radiative heating rate QJ

# F_RAD = np.zeros((lx, ly, lz))          # energy flux for RAD grid
F_int_x = np.zeros((lx, ly-1, lz-1))    # energy flux for inter face
F_int_y = np.zeros((lx-1, ly, lz-1))
F_int_z = np.zeros((lx-1, ly-1, lz))
Q_F = np.zeros((nx, ny, nz))            # radiative heating rate QF


### upwind point ###
# delta_l & dxc, dyc, dxd, dyd
for myu_num in range(0, 3):
    dl[myu_num] = min(D[0]/myu[myu_num][0], D[1]/myu[myu_num][1], D[2]/myu[myu_num][2])
    
    for dir in range(0, 3):
        dc[myu_num, dir] = myu[myu_num][dir] * dl[myu_num]
        dd[myu_num, dir] = D[dir] - dc[myu_num, dir]



##### define #####
# absorption coefficient for RAD grid 
def absorption_coefficient_RAD(alfa_MHD, i, j, k):
    if j == ly-1 and k == lz-1:
        return  exp((log(alfa_MHD[i,   j-1, k-1]) + log(alfa_MHD[i,   j-1, 0]) \
                   + log(alfa_MHD[i,   0,   k-1]) + log(alfa_MHD[i,   0,   0]) \
                   + log(alfa_MHD[i+1, j-1, k-1]) + log(alfa_MHD[i+1, j-1, 0]) \
                   + log(alfa_MHD[i+1, 0,   k-1]) + log(alfa_MHD[i+1, 0,   0])) / 8)              
            
    elif j == ly-1 and k != lz-1:
        return exp((log(alfa_MHD[i,   j-1, k-1]) + log(alfa_MHD[i,   j-1, k]) \
                  + log(alfa_MHD[i,   0,   k-1]) + log(alfa_MHD[i,   0,   k]) \
                  + log(alfa_MHD[i+1, j-1, k-1]) + log(alfa_MHD[i+1, j-1, k]) \
                  + log(alfa_MHD[i+1, 0,   k-1]) + log(alfa_MHD[i+1, 0,   k])) / 8)
                
    elif j != ly-1 and k == lz-1:
        return exp((log(alfa_MHD[i,   j-1, k-1]) + log(alfa_MHD[i,   j-1, 0]) \
                  + log(alfa_MHD[i,   j,   k-1]) + log(alfa_MHD[i,   j,   0]) \
                  + log(alfa_MHD[i+1, j-1, k-1]) + log(alfa_MHD[i+1, j-1, 0]) \
                  + log(alfa_MHD[i+1, j,   k-1]) + log(alfa_MHD[i+1, j,   0])) / 8)
            
    elif j == 0 and k == 0:
        return exp((log(alfa_MHD[i,   ny-1, nz-1]) + log(alfa_MHD[i,   ny-1, k]) \
                  + log(alfa_MHD[i,   j,    nz-1]) + log(alfa_MHD[i,   j,    k]) \
                  + log(alfa_MHD[i+1, ny-1, nz-1]) + log(alfa_MHD[i+1, ny-1, k]) \
                  + log(alfa_MHD[i+1, j,    nz-1]) + log(alfa_MHD[i+1, j,    k])) / 8)
            
    elif j == 0 and k != 0:
        return  exp((log(alfa_MHD[i,   ny-1, k-1]) + log(alfa_MHD[i,   ny-1, k]) \
                   + log(alfa_MHD[i,   j,    k-1]) + log(alfa_MHD[i,   j,    k]) \
                   + log(alfa_MHD[i+1, ny-1, k-1]) + log(alfa_MHD[i+1, ny-1, k]) \
                   + log(alfa_MHD[i+1, j,    k-1]) + log(alfa_MHD[i+1, j,    k])) / 8)
            
    elif j != 0 and k == 0:
        return exp((log(alfa_MHD[i,   j-1, nz-1]) + log(alfa_MHD[i,   j-1, k]) \
                  + log(alfa_MHD[i,   j,   nz-1]) + log(alfa_MHD[i,   j,   k]) \
                  + log(alfa_MHD[i+1, j-1, nz-1]) + log(alfa_MHD[i+1, j-1, k]) \
                  + log(alfa_MHD[i+1, j,   nz-1]) + log(alfa_MHD[i+1, j,   k])) / 8)

    else:
        return exp((log(alfa_MHD[i,   j-1, k-1]) + log(alfa_MHD[i,   j-1, k]) \
                  + log(alfa_MHD[i,   j,   k-1]) + log(alfa_MHD[i,   j,   k]) \
                  + log(alfa_MHD[i+1, j-1, k-1]) + log(alfa_MHD[i+1, j-1, k]) \
                  + log(alfa_MHD[i+1, j,   k-1]) + log(alfa_MHD[i+1, j,   k])) / 8)

# absorption coefficient for upwind point
def absorption_coefficient_up(alfa, i, j, k, dd, dc, D, myu_num):
    if j == 0 and k == 0:
        return exp((log(alfa[i,   j,    k   ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j,    k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   ly-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, ly-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j,    lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j,    lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i,   ly-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, ly-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])
                
    elif j != 0 and k == 0:
        return exp((log(alfa[i,   j,   k   ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j,   k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j,   lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j,   lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i,   j-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])

    elif j == 0 and k != 0:
        return exp((log(alfa[i,   j,    k  ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j,    k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   ly-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, ly-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j,    k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j,    k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i,   ly-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, ly-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])

    else:
        return exp((log(alfa[i,   j,   k  ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j,   k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i-1, j-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(alfa[i,   j,   k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j,   k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i,   j-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(alfa[i-1, j-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])


# delta opacity
def delta_opacity(alfa_up, i, j, k, dl, myu_num):

    if myu_num == 0:
        return alfa_up[myu_num, i, j, k] * (dl[myu_num]) * 4800000
    
    else:
        return alfa_up[myu_num, i, j, k] * (dl[myu_num]) * 4800000


# planck function for RAD grid 
def planck_function_RAD(B_MHD, i, j, k):
    if j == ly-1 and k == lz-1:
        return exp((log(B_MHD[i,   j-1, k-1]) + log(B_MHD[i,   j-1, 0]) \
                  + log(B_MHD[i,   0,   k-1]) + log(B_MHD[i,   0,   0]) \
                  + log(B_MHD[i+1, j-1, k-1]) + log(B_MHD[i+1, j-1, 0]) \
                  + log(B_MHD[i+1, 0,   k-1]) + log(B_MHD[i+1, 0,   0])) / 8)              
            
    elif j == ly-1 and k != lz-1:
        return exp((log(B_MHD[i,   j-1, k-1]) + log(B_MHD[i,   j-1, k]) \
                  + log(B_MHD[i,   0,   k-1]) + log(B_MHD[i,   0,   k]) \
                  + log(B_MHD[i+1, j-1, k-1]) + log(B_MHD[i+1, j-1, k]) \
                  + log(B_MHD[i+1, 0,   k-1]) + log(B_MHD[i+1, 0,   k])) / 8)
                
    elif j != ly-1 and k == lz-1:
        return exp((log(B_MHD[i,   j-1, k-1]) + log(B_MHD[i,   j-1, 0]) \
                  + log(B_MHD[i,   j,   k-1]) + log(B_MHD[i,   j,   0]) \
                  + log(B_MHD[i+1, j-1, k-1]) + log(B_MHD[i+1, j-1, 0]) \
                  + log(B_MHD[i+1, j,   k-1]) + log(B_MHD[i+1, j,   0])) / 8)
            
    elif j == 0 and k == 0:
        return exp((log(B_MHD[i,   ny-1, nz-1]) + log(B_MHD[i,   ny-1, k]) \
                  + log(B_MHD[i,   j,    nz-1]) + log(B_MHD[i,   j,    k]) \
                  + log(B_MHD[i+1, ny-1, nz-1]) + log(B_MHD[i+1, ny-1, k]) \
                  + log(B_MHD[i+1, j,    nz-1]) + log(B_MHD[i+1, j,    k])) / 8)
            
    elif j == 0 and k != 0:
        return exp((log(B_MHD[i,   ny-1, k-1]) + log(B_MHD[i,   ny-1, k]) \
                  + log(B_MHD[i,   j,    k-1]) + log(B_MHD[i,   j,    k]) \
                  + log(B_MHD[i+1, ny-1, k-1]) + log(B_MHD[i+1, ny-1, k]) \
                  + log(B_MHD[i+1, j,    k-1]) + log(B_MHD[i+1, j,    k])) / 8)
            
    elif j != 0 and k == 0:
        return exp((log(B_MHD[i,   j-1, nz-1]) + log(B_MHD[i,   j-1, k]) \
                  + log(B_MHD[i,   j,   nz-1]) + log(B_MHD[i,   j,   k]) \
                  + log(B_MHD[i+1, j-1, nz-1]) + log(B_MHD[i+1, j-1, k]) \
                  + log(B_MHD[i+1, j,   nz-1]) + log(B_MHD[i+1, j,   k])) / 8)

    else:
        return exp((log(B_MHD[i,   j-1, k-1]) + log(B_MHD[i,   j-1, k]) \
                  + log(B_MHD[i,   j,   k-1]) + log(B_MHD[i,   j,   k]) \
                  + log(B_MHD[i+1, j-1, k-1]) + log(B_MHD[i+1, j-1, k]) \
                  + log(B_MHD[i+1, j,   k-1]) + log(B_MHD[i+1, j,   k])) / 8)


# planck function for upwind point
def planck_function_up(B, i, j, k, dd, dc, D, myu_num):
    if j == 0 and k == 0:
        return exp((log(B[i,   j,    k   ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j,    k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   ly-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, ly-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j,    lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j,    lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i,   ly-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, ly-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])
                
    elif j != 0 and k == 0:
        return exp((log(B[i,   j,   k   ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j,   k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j,   lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j,   lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i,   j-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])

    elif j == 0 and k != 0:
        return exp((log(B[i,   j,    k  ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j,    k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   ly-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, ly-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j,    k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j,    k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i,   ly-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, ly-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])
    
    else:
        return exp((log(B[i,   j,   k  ]) * dd[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j,   k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i-1, j-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                  + log(B[i,   j,   k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j,   k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i,   j-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                  + log(B[i-1, j-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])


# delta intensity
def delta_intensity(B, delta_tu, i, j, k, myu_num):
    return B[i+1, j, k] * (1 - exp(- delta_tu[myu_num, i, j, k]))


# intensity for upwind point
def intensity_up(I, i, j, k, dd, dc, D, myu_num):
    if myu_num == 0:
        if j == 0 and k == 0:
            return exp((log(I[myu_num, i-1, j,    k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j,    lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
                
        elif j != 0 and k == 0:
            return exp((log(I[myu_num, i-1, j,   k   ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j,   lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])

        elif j == 0 and k != 0:
            return exp((log(I[myu_num, i-1, j,    k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j,    k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
    
        else:
            return exp((log(I[myu_num, i-1, j,   k  ]) * dc[myu_num, 0] * dd[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j,   k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
    
    elif myu_num == 1:
        if j == 0 and k == 0:
            return exp((log(I[myu_num, i,   ly-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i,   ly-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
                
        elif j != 0 and k == 0:
            return exp((log(I[myu_num, i,   j-1, k   ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k   ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i,   j-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])

        elif j == 0 and k != 0:
            return exp((log(I[myu_num, i,   ly-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i,   ly-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
    
        else:
            return exp((log(I[myu_num, i,   j-1, k  ]) * dd[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k  ]) * dc[myu_num, 0] * dc[myu_num, 1] * dd[myu_num, 2] \
                      + log(I[myu_num, i,   j-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                / D[0] * D[1] * D[2])
    
    else:
        if j == 0 and k == 0:
            return exp((log(I[myu_num, i,   j,    lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j,    lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i,   ly-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
                
        elif j != 0 and k == 0:
            return exp((log(I[myu_num, i,   j,   lz-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j,   lz-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i,   j-1, lz-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, lz-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])

        elif j == 0 and k != 0:
            return exp((log(I[myu_num, i,   j,    k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j,    k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i,   ly-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, ly-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])
    
        else:
            return exp((log(I[myu_num, i,   j,   k-1]) * dd[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j,   k-1]) * dc[myu_num, 0] * dd[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i,   j-1, k-1]) * dd[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2] \
                      + log(I[myu_num, i-1, j-1, k-1]) * dc[myu_num, 0] * dc[myu_num, 1] * dc[myu_num, 2])\
                    / D[0] * D[1] * D[2])


# mean intensiry for RAD grid
def mean_intensity_RAD(I):
    return (I[0, :, :, :] + I[1, :, :, :] + I[2, :, :, :]) / 3

# mean intensity for cell center
def mean_intensity_cellcenter(J_RAD, i, j, k):
    return (J_RAD[i,   j, k] + J_RAD[i,   j+1, k] + J_RAD[i,   j, k+1] + J_RAD[i,   j+1, k+1] \
          + J_RAD[i+1, j, k] + J_RAD[i+1, j+1, k] + J_RAD[i+1, j, k+1] + J_RAD[i+1, j+1, k+1]) \
          / 8

# radiative heating rate QJ
def radiative_heating_rate_QJ(J_cen, alfa_MHD, B_MHD, i, j, k):
    return 4 * pi * alfa_MHD[i, j, k] * (J_cen[i, j, k] - B_MHD[i, j, k])


# energy flux for RAD grid
def energy_flux_RAD(myu, I):
    return 4 / 3 * pi * (myu[0][0] * I[0, :, :, :] + myu[1][0] * I[1, :, :, :] + myu[2][0] * I[2, :, :, :])

# energy flux for inter face
def energy_flux_interface(F_RAD, i, j, k, direction):
    if direction == 0:
        return (F_RAD[i, j, k] + F_RAD[i, j+1, k] + F_RAD[i, j, k+1] + F_RAD[i, j+1, k+1]) / 4
    
    elif direction == 1:
        return (F_RAD[i, j, k] + F_RAD[i+1, j, k] + F_RAD[i, j, k+1] + F_RAD[i+1, j, k+1]) / 4
    
    else:
        return (F_RAD[i, j, k] + F_RAD[i, j+1, k] + F_RAD[i+1, j, k] + F_RAD[i+1, j+1, k]) / 4

# radiative heating rate QF
def radiative_heating_rate_QF(F_int, i, j, k):
    return - ((F_int[0, i+1, j,   k] - F_int[0, i, j, k]) \
            + (F_int[1, i,   j+1, k] - F_int[1, i, j, k]) \
            + (F_int[2, i,   j,   k+1] - F_int[2, i, j, k]))




##### opacity #####
# density
Ro0, tmp, tmp = np.meshgrid(ro0, d.p['y'], d.p['z'], indexing='ij')
ro = Ro0 + ro1   

# absorption_coefficient_MHD 
alfa_MHD = ro * op

# adxorption_coefficient_RAD
print('### alfa ###')

for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            
            alfa[i, j, k] = absorption_coefficient_RAD(alfa_MHD, i, j, k)
            
            # print(i, j, k, alfa[i, j, k])


# absorption coefficient for upwind point
print('### alfa for upwind point ###')

for myu_num in range(0, 3):
    for i in range(1, lx):
        for j in range(0, ly):
            for k in range(0, lz):
                alfa_up[myu_num, i-1, j, k] = absorption_coefficient_up(alfa, i, j, k, dd, dc, D, myu_num)
                
                # print(myu_num, i-1, j, k, alfa_up[myu_num, i-1, j, k])


# delta opacity
print('### delta tau ###')

for myu_num in range(0, 3):
    for i in range(lx-2, -1, -1):
        for j in range(ly-1, -1, -1):
            for k in range(lz-1, -1, -1):
                
                delta_tu[myu_num, i, j, k] = delta_opacity(alfa_up, i, j, k, dl, myu_num)

                print(1, i, j, k, delta_tu[1, i, j, k])



##### Intensity #####
# temperature
T0, tmp, tmp = np.meshgrid(t0, d.p['y'], d.p['z'], indexing = 'ij')
t = T0 + t1

# planck_function_MHD
sig = 5.67 * 10**(-5)   # Stefan-Btzmann constant [erg cm^-2 deg^-4 s^-1]
B_MHD = (sig * t**4) / pi

# planck function_RAD
print('### planck function ###')

for i in range(0, lx):
    for j in range(0, ly):
        for k in range(0, lz):
            B[i, j, k] = planck_function_RAD(B_MHD, i, j, k)
            
            # print(i, j, k, B[i, j, k])


# planck_function_upwind
print('### planck function for upwind point ###')

for myu_num in range(0, 3):
    for i in range(1, lx):
        for j in range(0, ly):
            for k in range(0, lz):
                
                B_up[myu_num, i-1, j, k] = planck_function_up(B, i, j, k, dd, dc, D, myu_num)
                
                # print(myu_num, i, j, k, B_up[myu_num, i-1, j, k])


##### delta intensity #####
print('### delta intensity ###')

# x_bottom
lxb = 50

for myu_num in range(0, 3):
    for i in range(lxb, lx-1):
        for j in range(0, ly):
            for k in range(0, lz):
                
                delta_I[myu_num, i, j, k] = delta_intensity(B, delta_tu, i, j, k, myu_num)

                # print(delta_tu[myu_num, i, j, k])


##### intensity #####
print('### intensity ###')

# boundary condition
I[:, lxb, :, :] = B[lxb, :, :]
I[1,   :, ly-1, :] = B[:,   ly-1, :]
I[2,   :, :, lz-1] = B[:,   :, lz-1]

for myu_num in range(0, 3):
    # myu1
    if myu_num == 0:
        for i in range(lxb+1, lx):
            for j in range(0, ly):
                for k in range(0, lz):
                    # upwind point
                    I_up[myu_num, i-1, j, k] = intensity_up(I, i, j, k, dd, dc, D, myu_num)
                    I[myu_num, i, j, k] = I_up[myu_num, i-1, j, k] * exp(- delta_tu[myu_num, i-1, j, k]) + delta_I[myu_num, i-1, j, k]
                    # print(i, j, k, I[myu_num, i, j, k])
    
    # myu2
    elif myu_num == 1:
        while True:
            for j in range(0, ly):
                for i in range(lxb+1, lx):
                    for k in range(0, lz):
                        # upwind point
                        I_up[myu_num, i-1, j, k] = intensity_up(I, i, j, k, dd, dc, D, myu_num)
                        I[myu_num, i, j, k] = I_up[myu_num, i-1, j, k] * np.exp(-delta_tu[myu_num, i-1, j, k]) + delta_I[myu_num, i-1, j, k]
                        # print(i, j, k, I[myu_num, i, j, k])

            # error calculation and updating I_up
            max_error = 0
            for i in range(lxb+1, lx):
                for k in range(0, lz):
                    I_up_new[myu_num, i-1, 0, k] = intensity_up(I, i, 0, k, dd, dc, D, myu_num)
                    err = I_up[myu_num, i-1, 0, k] - I_up_new[myu_num, i-1, 0, k]

                    max_error = max(max_error, abs(err))

                    # Update I_up with I_up_new
                    I_up[myu_num, i-1, 0, k] = I_up_new[myu_num, i-1, 0, k]

            # Check if the maximum error is below the threshold
            if max_error <= 10**(-3):
                for j in range(0, ly):
                    for i in range(lxb+1, lx):
                        for k in range(0, lz):
                            # upwind point
                            I_up[myu_num, i-1, j, k] = intensity_up(I, i, j, k, dd, dc, D, myu_num)
                            I[myu_num, i, j, k] = I_up[myu_num, i-1, j, k] * np.exp(-delta_tu[myu_num, i-1, j, k]) + delta_I[myu_num, i-1, j, k]
                            # print(i, j, k, I[myu_num, i, j, k])
                break
    
    # myu3
    elif myu_num == 2:
        while True:
            for k in range(0, lz):
                for i in range(lxb+1, lx):
                    for j in range(0, ly):
                        # upwind point
                        I_up[myu_num, i-1, j, k] = intensity_up(I, i, j, k, dd, dc, D, myu_num)
                        I[myu_num, i, j, k] = I_up[myu_num, i-1, j, k] * np.exp(-delta_tu[myu_num, i-1, j, k]) + delta_I[myu_num, i-1, j, k]
                        # print(i, j, k, I[myu_num, i, j, k])

            # error calculation and updating I_up
            max_error = 0
            for i in range(lxb+1, lx):
                for j in range(0, ly):
                    I_up_new[myu_num, i-1, j, 0] = intensity_up(I, i, j, 0, dd, dc, D, myu_num)
                    err = I_up[myu_num, i-1, j, 0] - I_up_new[myu_num, i-1, j, 0]

                    max_error = max(max_error, abs(err))

                    # Update I_up with I_up_new
                    I_up[myu_num, i-1, j, 0] = I_up_new[myu_num, i-1, j, 0]

            # Check if the maximum error is below the threshold
            if max_error <= 10**(-3):
                for k in range(0, lz):
                    for i in range(lxb+1, lx):
                        for j in range(0, ly):
                            # upwind point
                            I_up[myu_num, i-1, j, k] = intensity_up(I, i, j, k, dd, dc, D, myu_num)
                            I[myu_num, i, j, k] = I_up[myu_num, i-1, j, k] * np.exp(-delta_tu[myu_num, i-1, j, k]) + delta_I[myu_num, i-1, j, k]
                            # print(i, j, k, I[myu_num, i, j, k])
                break


##### radiative heating rate #####
### radiative heating rate QJ ###
print('### radiative heating rate QJ ###')

# x_upper
nxu = 52     # x_upper is point of tu < 0.1

# mean intensity for RAD grid
J_RAD = mean_intensity_RAD(I)
# print(J_RAD)

# mean intensity for cell center
for i in range(lxb, nxu+1):
    for j in range(0, ny):
        for k in range(0, nz):
            J_cen[i, j, k] = mean_intensity_cellcenter(J_RAD, i, j, k)
            # print(i, j, k, J_cen[i, j, k])

# radiative heating rate QJ 
for i in range(lxb, nxu+1):
    for j in range(0, ny):
        for k in range(0, nz):
            Q_J[i, j, k] = radiative_heating_rate_QJ(J_cen, alfa_MHD, B_MHD, i, j, k)
            # print(i, j, k, Q_J[i, j, k])


### radiative heating rate QF ###
print('### radiative heating rate QF ###')

# energy flux for RAD grid
F_RAD = energy_flux_RAD(myu, I)
# print(F_RAD)

# energy flux for inter face
for direction in range(0, 3):
    if direction == 0:
        for i in range(nxu, lx):
            for j in range(0, ly-1):
                for k in range(0, lz-1):
                    F_int_x[i, j, k] = energy_flux_interface(F_RAD, i, j, k, direction)
                    # print(direction, F_int_x[i, j, k])
    
    elif direction == 1:
        for i in range(nxu, lx-1):
            for j in range(0, ly):
                for k in range(0, lz-1):
                    F_int_y[i, j, k] = energy_flux_interface(F_RAD, i, j, k, direction)
                    print(direction, F_int_y[i, j, k])
    
    else:
        for i in range(nxu, lx-1):
            for j in range(0, ly-1):
                for k in range(0, lz):
                    F_int_z[i, j, k] = energy_flux_interface(F_RAD, i, j, k, direction)
                    print(direction, F_int_z[i, j, k])

# radiative heating rate QF 
for i in range(nxu, nx-2):
    for j in range(0, ny):
        for k in range(0, nz):
            Q_F[i, j, k] = radiative_heating_rate_QF(F_int, i, j, k)
            print(i, j, k, Q_F[i, j, k])


##### plot #####
print('### plot ###')
Intensity = I[1, 62, :, :]
# ymax = (ly-1) * myu[2][1]
# zmax = (lz-1) * myu[2][2]

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure()
plt.imshow(Intensity)   # extent = (0, ymax, 0, zmax)

ticks = np.linspace(Intensity.min(), Intensity.max(), 10)
cbar = plt.colorbar()
cbar.set_ticks(ticks)

save_dir = '/cidashome/sc/c0287shimizu/R2D2_ex/py/fig/MR/MR_test'
file_name = 'test_old.png'
file_path = os.path.join(save_dir, file_name)

plt.savefig(file_path)
plt.close()