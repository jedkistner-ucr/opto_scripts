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

path = "C:/Jed/Built"
name = "2021_11_03_2_sdcube" #WS2-WSe2 sdcube
# name = "2023_03_29_sdcube" #ao07 sd cube
# name = "2022_02_04_11_sdcube" #d02 sd cube
# name = '2022_01_27_0_sdcube'   #sd cube of ao02 at -5.25
# name = '2022_01_27_1_sdcube'  #sd cube of ao02 at -4.48
# name = "2022_02_03_33_powercube" #tg -bg powercube for Ao02
# name = "2022_01_26_65_powercube" #tg -bg powercube for Ao02
# name = "2022_02_03_12_powercube"  #sd_gate powercube ao07
# name = "2022_02_03_12_powercube_butter_5_2"  #sd-gate powercube for ao07
# name = '2021_03_25_51_powercube'
# name = '2021_03_25_51_powercube_butter_4_2.0'
xval, yval, cval, data, rf, powdata, info_ = loader.load_cube(path, name)

# name = "2022_01_27_powermap_noa"
# xval, yval, cval, zval, data, rf, power = loader.load_hypercube(path, name)
# s = 25
# data = data[:,:,:,s]
# print(zval[s])

# data = (data*1e9)+2.8 #for 2022_01_27_0
# data = rf
# data = (data*1e9)+.8
# data = (data*1e9)+20.8 #for 2022_02_03_12
data = (data*1e9)

use_power_for_cval = False
cubelabel = "Vsd (V)"
# use_power_for_cval = True
# cubelabel = "Power (mW)"

save = True
show = True
savename = name + "0"
normalize = True
# normalize = True
contours = True
contours = False
rfoverlay = False
zmin = -10
zmax = 10
# zmin = 0
# zmax = 10
zmin = 0
zmax = 0

conduction = True
secondder = False
symmetric_scale = True
scale = .9
smooth = 0

mapcolor = 'seismic'
# mapcolor = 'plasma'
# mapcolor = "Reds"
wavelength = 810

spotsmooth = 0
spot = [-16.73, -39.12] #black
spot2 = [-17, -45] #blue
spot3 = [] #red
spot = [17.7, -13.5]
spot = [22.9, -11.5] #black
spot = [19.8, -10] #black
# spot2 = [28.8, -11.3] #blue
# spot3 = [23, -5.3] #red
# spot = [.85, -5.5]
# spot2 = [.3, -5]
# spot3 = [.4, -1.4]
# spot = [-6.44, -4]
# spot2 = [-3.74, -6.57]
# spot3 = [-2.25, -2.25]
spot = []
spot2 = []
spot3 = []

vline = [1, -1, 0.22, -.2]
hline = []

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
if contours:
        savename = savename + "_rfcont"
if rfoverlay:
        savename = savename + "_rfover"

# print(cval.shape)
# print(zval.shape)
clen = cval.size
xlen = xval.size
ylen = yval.size

xlow = 0
xhigh = 0
ylow = 0
yhigh = 0

if use_power_for_cval:
        cval = []
        for i in range(clen):
                cval.append(np.mean(powdata[:,:,i]))
        cval = np.asarray(cval)
        print('Bias values: ')
        print(cval)


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


if xval[0] > xval[1] and yval[0] > yval[1]:
        xt = get_extent(xval, yval,invertx= True ,inverty= False)
elif xval[0] > xval[1] and yval[0] < yval[1]:
        xt = get_extent(xval, yval,invertx= True ,inverty= True)
elif xval[0] < xval[1] and yval[0] > yval[1]:
        xt = get_extent(xval, yval,invertx= False ,inverty= False)
else:
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

# spot0_ = np.searchsorted(-xval, -spot[0])
# spot1_ = np.searchsorted(-yval, -spot[1])
# for i in range(clen):
#         data[spot1_, spot0_,i] = np.max(data)*2


for i in range(clen):

                d = data[:,:,i]
                rfd = rf[:,:,i]

                if smooth > 0:
                        d = applyGauss(d, smooth)

                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)
                if rfoverlay:
                        rcmap, rcnorm = color_memap("gray", rfd)
                        ax[i].imshow(rfd, cmap = rcmap, norm = rcnorm, extent = xt, aspect = 'auto', origin = 'lower', alpha = 1)
                if rfoverlay:
                        ax[i].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower', alpha = 0.4)
                else:
                        ax[i].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[i].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 4, aspect = 'auto', origin = 'lower')

                # if p == sdindex0:
                #         ax[axi][axp].set_ylabel(str(np.around(zval[i], 2)) + "V g")
                ax[i].text(.06, .95, str(np.around( cval[i], 3) ) + cubelabel , horizontalalignment = 'left', verticalalignment = 'top', transform=ax[i].transAxes)

                for v in vline:
                        ax[i].axvline(v, c = 'k', linestyle = ":")

fig.suptitle(name + "  ||  PC images with varying value...")

if len(spot) > 0:
        axs, figs = makeaxes(1)
        fig1 = pl.figure(figsize=figs)
        ax1 = fig1.add_axes(axs[0])

        if xval[0] < xval[1]:
                spot[0] = np.searchsorted(xval, spot[0])
        else:
                spot[0] = np.searchsorted(-xval, -spot[0])
        if yval[0] < yval[1]:
                spot[1] = np.searchsorted(yval, spot[1])
        else:
                spot[1] = np.searchsorted(-yval, -spot[1])

        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot[1]-1+m, spot[0]-1-n,:]

        # y += data[spot[1], spot[0],:]
        y = y / 9
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'k')
        for i in range(clen):
                ax[i].scatter(xval[spot[0]], yval[spot[1]], color = 'k', alpha = .5)
        ax1.set_xlim(cval[0], cval[-1])
        ax1.set_xlabel(cubelabel)
        ax1.set_ylabel("Current (nA)")

if len(spot2) > 0:        
        if xval[0] < xval[1]:
                spot2[0] = np.searchsorted(xval, spot2[0])
        else:
                spot2[0] = np.searchsorted(-xval, -spot2[0])
        if yval[0] < yval[1]:
                spot2[1] = np.searchsorted(yval, spot2[1])
        else:
                spot2[1] = np.searchsorted(-yval, -spot2[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot2[1]-1+m, spot2[0]-1-n,:]
        y = y / 9
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'b')
        for i in range(clen):
                ax[i].scatter(xval[spot2[0]], yval[spot2[1]], color = 'b', alpha = .5)

if len(spot3) > 0:        
        if xval[0] < xval[1]:
                spot3[0] = np.searchsorted(xval, spot3[0])
        else:
                spot3[0] = np.searchsorted(-xval, -spot3[0])
        if yval[0] < yval[1]:
                spot3[1] = np.searchsorted(yval, spot3[1])
        else:
                spot3[1] = np.searchsorted(-yval, -spot3[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += data[spot3[1]-1+m, spot3[0]-1-n,:]
        y = y / 9
        if spotsmooth > 0:
                y = ndimage.gaussian_filter1d(y, spotsmooth)
        ax1.plot(x, y, linewidth = 3, c = 'r')
        for i in range(clen):
                ax[i].scatter(xval[spot3[0]], yval[spot3[1]], color = 'r', alpha = .5)

savename = savename + "_cube_"

print('saving')

if save:
        fig.savefig("C:/Jed/Figures/" + savename + ".png")
        if len(spot) > 0:
                fig1.savefig("C:/Jed/Figures/" + savename + 'fig1' + ".png")
if show:
        pl.show()