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
name = "2023_03_28_powermap_noa"
name = "2022_01_26_biasmap_noa" 
xval, yval, cval, zval, data, rf, power = loader.load_hypercube(path, name)
data = (data*1e9)

display_bias = -1
usepowerforcval = False

save = False
show = True
savename = name + ""
normalize = True
contours = False
rfoverlay = False
zmin = 0
zmax = 0

secondder = False
smooth = 0

mapcolor = 'seismic'
mapcolor = 'plasma'
wavelength = 810

spot = [-16.73, -39.12]
spot2 = [-17, -45]
spot = [22.1, -3]
spot2 = []

# gate = [0, 2]
# source = [3, 3]


# if conduction:
#         derstep = np.abs(cval[0] - cval[1])
#         data = np.gradient(data, derstep, axis = 2)
#         if secondder:
#                 data = np.gradient(data, derstep, axis = 2)
#                 savename = savename + "secondDerivative"
#         else:
#                 savename = savename + "conduction"


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

if usepowerforcval:
        cval = []
        for i in range(clen):
                cval.append(np.mean(power[:,:,i,0]))
        cval = np.asarray(cval)
print('Bias values: ')
print(zval)

biasdisplayindex = np.searchsorted(zval, display_bias)
if biasdisplayindex == zlen:
        biasdisplayindex -= 1

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

subdata = data[:,:,:, biasdisplayindex]
subrf = rf[:,:,:,  biasdisplayindex]

xt = get_extent(xval, yval, inverty= False)

mg = 20
if normalize:
        
        zmax = np.max(subdata)
        zmin = np.min(subdata)

print('populating')

for i in range(clen):

                d = subdata[:,:,i]
                rfd = subrf[:,:,i]

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
                # ax[axi][axp].text(.06, .95, str(np.around( cval[p], 3) ) + "V sd" , horizontalalignment = 'left', verticalalignment = 'top', transform=ax[axi][axp].transAxes)

fig.suptitle(name + "  ||  PC images with varying bias at " + str(wavelength) + "nm and" + str(display_bias) + "(V)" )

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
                        y += subdata[spot[1]-1+m, spot[0]-1-n,:]
        ax1.plot(x, y, linewidth = 2, c = 'k')
        for i in range(clen):
                ax[i].scatter(xval[spot[0] -1 + m], yval[spot[1] -1 + n], color = 'k', alpha = .5)
        ax1.set_xlim(cval[0], cval[-1])
        ax1.set_xlabel("Power (mW)")
        ax1.set_ylabel("Current (nA)")

if len(spot2) > 0:        
        spot2[0] = np.searchsorted(xval, spot2[0])
        spot2[1] = np.searchsorted(yval, spot2[1])
        x = cval
        y = np.zeros(cval.shape)
        for m in range(3):
                for n in range(3):
                        y += subdata[spot2[1]-1+m, spot2[0]-1-n,:]
        ax1.plot(x, y, linewidth = 2, c = 'b')
        for i in range(clen):
                ax[i].scatter(xval[spot2[0] -1 + m], yval[spot2[1] -1 + n], color = 'b', alpha = .5)

savename = savename + "_subset_"

print('saving')

if save:
        fig.savefig("C:/Jed/Figures/" + savename + ".png")
if show:
        pl.show()