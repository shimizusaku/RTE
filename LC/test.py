# This code is RTE for multi rays.
import os
import time
# import R2D2
import numpy as np
# import matplotlib.pyplot as plt
from math import log, exp, pi, copysign

# start_time = time.time()

# ##### parameter #####
# ### R2D2 data ###
# d = R2D2.R2D2_data('/mnt/solar07a/c0287shimizu/R2D2/run/d001/data/')
# d.read_qq(10,'ro')    
# d.read_qq(10,'op')    
# d.read_qq(10, 'te')  
# d.read_qq(10, 'se')  

# ro0 = d.p['ro0']      # background field density
# ro1 = d.qq['ro']      # density disturbances
# op = d.qq['op']       # opacity
# x = d.p['x']          # height of index i
# y = d.p['y']          # lenght of index j
# z = d.p['z']          # lenght of index k
# t1 = d.qq['te']       # temperature disturbances
# t0 = d.p['te0']       # background of tenperature
# se0 = d.p['se0']      # background field entropy
# se1 = d.qq['se']      # entropy disturbances


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



print('dxc:', dxc[:])
print('dyc:', dyc[:])
print('dzc:', dzc[:])