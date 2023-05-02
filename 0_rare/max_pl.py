'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

# import matplotlib as mpl
# import matplotlib.pyplot as pl
# import matplotlib.colors as colors
# from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

import sys
# sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map


name = "d01_m9_5K_plstack_tg_bg"
name = "d01_m9_5K_plstack_tg_neg_bg"
# name = "d02_deltav"
# name = "d02_sd"
# name = "d02_totalv"
# name = "d02_deltav"
# name = "A02_2gate_sd"
# name = "A02_0gate_widesd"
# name = 'ao2_totalv'
name = 'ao2_deltaV'
# name = 'ao3_deltav'
# name = 'ao4_deltav'
# name = 'ao1_deltaV'
# name = 'ao1_totalV'

poweraxis = ""

path = "C:/Jed/Built"
xval, yval, data = loader.load_simplemap(path, name)

xlabel = ""
ylabel = ''

offset = 0
# data = np.flip(data, 1)
# data = data * -1

save = True
show = True
extraplot = False

#Parameters
mapcolor_name = 'plasma' #Color of map
mapcolor_name = 'magma' #Color of map
mapcolor_name = "viridis"
# mapcolor_name = 'seismic'
colormap_name = 'viridis_r' #Color of slices

ylen = yval.size
xlen = xval.size

#Loads and builds axes
axs, figs = makeaxes(2)
ax = []
axsc = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

toE = ((4.4/3.4)*(30))+1.2

maxval = np.zeros(yval.shape)
for i in range(ylen):
    maxval[i] = np.max(data[i,:])

# yval = yval/toE
ax[1].scatter(yval, maxval, color = 'k')


# if save:
#     pl.savefig("C:/Jed/Figures/" + savename + ".png")

if show:
    pl.show()