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
# name = "2020_08_12_biascube"
name = "2020_09_29_biascube_noa"
# name = "2020_09_03_biascube_noa"
name = "2020_10_27_noa_overlap"
name = "2021_03_25_biasmap_noa"
name = "2021_03_24_powermap_noa"
# name = "2023_03_28_powermap_noa"
xval, yval, cval, zval, data, rf, power = loader.load_hypercube(path, name)
data = (data*1e12)
# cval = np.linspace(.7, 2.3, 50)

# zval = np.flip(zval, 0)
# data = np.flip(data, 3)
# cval = cval * -1

save = False
show = True
savename = name + ""
normalize = True
normalizegate = False
contours = True
rfoverlay = False
zmin = 0
zmax = 0
# zmin = -.1
# zmax = .1
# zmin = -.03
# zmax = .017

conduction = False
secondder = False
smooth = False
smoothval = 1

mapcolor = 'seismic'
mapcolor = 'plasma'
wavelength = 850

gate = [0, 2]
source = [3, 3]

if normalizegate:
        normalize = False

if conduction:
        derstep = np.abs(cval[0] - cval[1])
        data = np.gradient(data, derstep, axis = 2)
        if secondder:
                data = np.gradient(data, derstep, axis = 2)
                savename = savename + "secondDerivative"
        else:
                savename = savename + "conduction"

if normalizegate:
        savename = savename + "_normgate"
if normalize:
        savename = savename + "_norm"
if contours:
        savename = savename + "_rfcont"
if rfoverlay:
        savename = savename + "_rfover"

# print(cval.shape)
# print(zval.shape)
clen = cval[:].size
zlen = zval[:].size
xlen = xval[:].size
ylen = yval[:].size

xlow = 0
xhigh = 0
ylow = 0
yhigh = 0

gateindex0 = np.searchsorted(zval, gate[0])
gateindex1 = np.searchsorted(zval, gate[1])
sdindex0 = np.searchsorted(cval, source[0])
sdindex1 = np.searchsorted(cval, source[1])
sdindex0 = -2
sdindex1 = -1

if gateindex0 == zlen:
        gateindex0 = gateindex0 - 1
if gateindex1 == zlen:
        gateindex1 = gateindex1 - 1
if sdindex0 == clen:
        sdindex0 = sdindex0 - 1
if sdindex1 == clen:
        sdindex1 = sdindex1 - 1

ncols = np.abs(sdindex0 - sdindex1 )+1
nrows = np.abs(gateindex0 - gateindex1)+1

print('generating plots')

axs, figs = makemapgrid(nrows,ncols, scale = 1)
fig = pl.figure(figsize=figs)
ax = []
for h in axs:
        k = []
        for a in h:
                k.append(fig.add_axes(a))
        ax.append(k)

for k in ax:
        for a in k:
                a.set_xticks([])
                a.set_yticks([])

xt = get_extent(xval, yval, inverty= False)

mg = 20
if normalize:
        
        zmax = np.max(data[mg:ylen-mg,mg:xlen-mg,sdindex0:sdindex1,gateindex0:gateindex1])
        zmin = np.min(data[mg:ylen-mg,mg:xlen-mg,sdindex0:sdindex1,gateindex0:gateindex1])

print('populating')

for i in range(gateindex0, gateindex1+1):
        if normalizegate:
                # zmax = np.max(data[mg:ylen-mg,mg:xlen-mg,sdindex0:sdindex1,i])
                # zmin = np.min(data[mg:ylen-mg,mg:xlen-mg,sdindex0:sdindex1,i])
                zmax = np.max(data[mg:ylen-mg,mg:xlen-mg,:,i])
                zmin = np.min(data[mg:ylen-mg,mg:xlen-mg,:,i])
        for p in range(sdindex0, sdindex1):

                axi = i - gateindex0
                axp = p - sdindex0

                d = data[:,:,p,i]
                rfd = rf[:,:,p,i]

                if smooth:
                        d = applyGauss(d, smoothval)

                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)
                if rfoverlay:
                        rcmap, rcnorm = color_memap("gray", rf[mg:ylen-mg,mg:xlen-mg,p,i])
                        ax[axi][axp].imshow(rf[:,:,p,i], cmap = rcmap, norm = rcnorm, extent = xt, aspect = 'auto', origin = 'lower', alpha = 1)
                if rfoverlay:
                        ax[axi][axp].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower', alpha = 0.4)
                else:
                        ax[axi][axp].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[axi][axp].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 4, aspect = 'auto', origin = 'lower')

                if p == sdindex0:
                        ax[axi][axp].set_ylabel(str(np.around(zval[i], 2)) + "V g")
                ax[axi][axp].text(.06, .95, str(np.around( cval[p], 3) ) + "V sd" , horizontalalignment = 'left', verticalalignment = 'top', transform=ax[axi][axp].transAxes)

fig.suptitle(name + "  ||  PC images with varying bias at " + str(wavelength) + "nm")

savename = savename + "_subset_"+str(zval[gateindex0]) +str(zval[gateindex1])+str(cval[sdindex0])+str(cval[sdindex1])

print('saving')

if save:
        pl.savefig("C:/QMO/generated plots/" + savename + ".png")
if show:
        pl.show()