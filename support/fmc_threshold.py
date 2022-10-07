# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:18:36 2022

@author: mc16535
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import scipy.interpolate as sit
import scipy.io as sio
import os

os.chdir(r"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE")
mat = sio.loadmat("a_0mm_d_1_b_BP.mat")

data = np.abs(mat['data']).T
time = mat['time'][:, 0]

# Reduce resolution
grid_x, grid_y   = np.meshgrid(time, range(data.shape[0]))
ngrid_x, ngrid_y = np.meshgrid(time[::2], range(data.shape[0])[::4])
pts = np.asarray([grid_x.ravel(), grid_y.ravel()]).T
vls = data.ravel()
grid_z = sit.griddata(pts, vls, (ngrid_x, ngrid_y), method='nearest')

max_, min_ = np.max(grid_z)/10, np.min(grid_z)
threshold  = (max_ + min_) / 2
mask = np.ones(grid_z.shape)
mask[grid_z < threshold] = 0

fig = plt.figure()
plt.subplots_adjust(bottom=0.25)
ax = fig.subplots()
p = ax.imshow(mask, cmap="gist_yarg", extent=[time[0], time[-1], data.shape[1], 1], aspect='auto')

ax_slide = plt.axes([0.25, 0.1, 0.65, 0.03])

frame_slider = Slider(ax_slide, 'threshold',
                  min_, max_, valinit=(max_ + min_) / 2, valstep=(max_ + min_) / 50)

def update_frame(val):
    threshold = frame_slider.val
    mask[grid_z >= threshold] = 1
    mask[grid_z < threshold]  = 0
    p.set_data(mask)
    fig.canvas.draw()
    
frame_slider.on_changed(update_frame)
plt.show()



# fig = plt.figure()
# plt.subplots_adjust(bottom=0.25)
# ax = fig.subplots()
# p = ax.imshow(data[0, :, :], extent=[xlims[0], xlims[1], ylims[0], ylims[1]], vmin=0, vmax=1e-11)
 
# ax_slide = plt.axes([0.25, 0.1, 0.65, 0.03])
 
# frame_slider = Slider(ax_slide, 'Frame',
#                   0, 160, valinit=0, valstep=1)

# def update_frame(val):
#     current_v = frame_slider.val
#     p.set_data(data[current_v, :, :])
#     fig.canvas.draw()
 
frame_slider.on_changed(update_frame)
plt.show()