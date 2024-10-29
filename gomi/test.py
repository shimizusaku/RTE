import os
import time
import R2D2
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, pi

start_time = time.time()

##### parameter #####
### R2D2 data ###
d = R2D2.R2D2_data('../../run/d007/data/')

d.read_qq_rt(5, 'in')
Intensity = d.qd['in']
qrad = d.qd['qrad']

x_point = [1, 62]
mti = ['0', '1' , '2']
iti = ['0', '1']
jti = ['0', '1']
kti = ['0', '1']

for m in range(0, 3):
    for idir in range(0, 2):
        for jdir in range(0, 2):
            for kdir in range(0, 2):
                plot = Intensity[m, idir, jdir, kdir, x_point[idir], :, :]

                plt.rcParams['xtick.direction'] = 'in'
                plt.rcParams['ytick.direction'] = 'in'

                plt.figure()
                plt.imshow(plot)
                plt.title(mti[m] + iti[idir] + jti[jdir] + kti[kdir])

                ticks = np.linspace(plot.min(), plot.max(), 10)
                cbar = plt.colorbar()
                cbar.set_ticks(ticks)

                save_dir = '/cidashome/sc/c0287shimizu/R2D2_ex/py/fig/MR/d007/hotta_RTE'
                file_name = mti[m] + iti[idir] + jti[jdir] + kti[kdir] + '.png'
                file_path = os.path.join(save_dir, file_name)

                plt.savefig(file_path)
                plt.close


x_point_ = [1, 46, 62]
title = ['bottom', 'tu1', 'top']

for i in range(0, 3):
    plot = qrad[x_point_[i], :, :]

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.figure()
    plt.imshow(plot)
    plt.title('Qrad ' + title[i])

    ticks = np.linspace(plot.min(), plot.max(), 10)
    cbar = plt.colorbar()
    cbar.set_ticks(ticks)

    save_dir = '/cidashome/sc/c0287shimizu/R2D2_ex/py/fig/MR/d007/hotta_RTE'
    file_name = 'Qrad_' + title[i] + '.png'
    file_path = os.path.join(save_dir, file_name)

    plt.savefig(file_path)
    plt.close