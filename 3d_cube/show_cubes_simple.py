'''
Takes a built bias hypecube and works with it
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp

from datetime import date

from plot_tools import *
from analysis_tools import *
import loader
# from builder import build_map

path = "C:/QMO/Built"

name = 'erfu01_sd_bg'
# name = 'erfu01_tg_bg'
xval0, yval0, cval0, data = loader.load_simplecube(path, name)
#for erfu's maps cval is wavelength and yval is gate (sd is either the other gate or source drain)
#i just want to rearrange these a little
xval, yval, cval = cval0, xval0, yval0
data = np.swapaxes(data, 0, 2)
data = np.swapaxes(data, 2, 1)
data = np.swapaxes(data, 0, 1)


use_power_for_cval = False
cubelabel = "Vbg (V)"
# use_power_for_cval = True
# cubelabel = "Power (mW)"

save = True
show = True
savename = name + ""
normalize = False
# normalize = True
zmin = -10
zmax = 10
zmin = 0
zmax = 0

conduction = False
secondder = False
symmetric_scale = False
scale = .9
smooth = 0

mapcolor = 'seismic'
# mapcolor = 'plasma'
mapcolor = "Reds"
wavelength = 810

spotsmooth = 0
spot = [-16.73, -39.12] #black
spot2 = [-17, -45] #blue
spot3 = [] #red
spot = [17.7, -13.5]
spot = [22.9, -11.5] #black
spot2 = [28.8, -11.3] #blue
spot3 = [23, -5.3] #red
spot = []
spot2 = []
spot3 = []

if conduction:
        derstep = np.abs(cval[0] - cval[1])
        data = np.gradient(data, derstep, axis = 2)
        if secondder:
                data = np.gradient(data, derstep, axis = 2)
                savename = savename + "secondDerivative"
        else:
                savename = savename + "conduction"

if normalize:
        savename = savename + "_norm"

# print(cval.shape)
# print(zval.shape)
clen = cval.size
xlen = xval.size
ylen = yval.size

xlow = 0
xhigh = 0
ylow = 0
yhigh = 0

print('generating plots')

if clen % 10 == 0:
        nrows = clen // 10
else:
        nrows = (clen//10)+1

axs, figs = makemapgrid(nrows,10, scale = 1)
fig = pl.figure(figsize=figs)
ax = []
for h in axs:
        for a in h:
                ax.append(fig.add_axes(a))

for a in ax:
        a.set_xticks([])
        a.set_yticks([])

xt = get_extent(xval, yval, inverty= True)

mg = 20
if normalize:
        zmax = np.max(data)
        zmin = np.min(data)

if symmetric_scale:
        mapcolor = 'bwr'
        zmax = np.max(data) * scale
        zmin = -zmax * scale

print('populating')

# spot0_ = np.searchsorted(xval, spot[0])
# spot1_ = np.searchsorted(yval, spot[1])
# for i in range(clen):
#         data[spot1_, spot0_,i] = np.max(data)*2


for i in range(clen):

                d = data[:,:,i]

                if smooth > 0:
                        d = applyGauss(d, smooth)

                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)
                ax[i].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')

                ax[i].text(.06, .95, str(np.around( cval[i], 3) ) + cubelabel , horizontalalignment = 'left', verticalalignment = 'top', transform=ax[i].transAxes)

fig.suptitle(name + "  ||  images with varying value...")

if len(spot) > 0:
        axs, figs = makeaxes(1)
        fig1 = pl.figure(figsize=figs)
        ax1 = fig1.add_axes(axs[0])

        spot[0] = np.searchsorted(xval, spot[0])
        spot[1] = np.searchsorted(yval, spot[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot[1]-1+m, spot[0]-1-n,:]
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'k')
        for i in range(clen):
                ax[i].scatter(xval[spot[0] -1 + m], yval[spot[1] -1 + n], color = 'k', alpha = .5)
        ax1.set_xlim(cval[0], cval[-1])
        ax1.set_xlabel(cubelabel)
        ax1.set_ylabel("Current (nA)")

if len(spot2) > 0:        
        spot2[0] = np.searchsorted(xval, spot2[0])
        spot2[1] = np.searchsorted(yval, spot2[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot2[1]-1+m, spot2[0]-1-n,:]
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'b')
        for i in range(clen):
                ax[i].scatter(xval[spot2[0] -1 + m], yval[spot2[1] -1 + n], color = 'b', alpha = .5)

if len(spot3) > 0:        
        spot3[0] = np.searchsorted(xval, spot3[0])
        spot3[1] = np.searchsorted(yval, spot3[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot3[1]-1+m, spot3[0]-1-n,:]
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'r')
        for i in range(clen):
                ax[i].scatter(xval[spot3[0] -1 + m], yval[spot3[1] -1 + n], color = 'r', alpha = .5)

savename = savename + "_cube_"

print('saving')

if save:
        fig.savefig("C:/QMO/generated plots/" + savename + ".png")
        if len(spot) > 0:
                fig1.savefig("C:/QMO/generated plots/" + savename + 'fig1' + ".png")
if show:
        pl.show()