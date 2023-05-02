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
name = "2021_10_29_powermap_noa"  
name = "2022_01_29_powermap"  
name = "2022_01_26_biasmap_noa" 
name = "2021_03_25_biasmap_noa"
name = "2023_01_30_biasmap_noa"  
# name = "2021_03_24_powermap_noa"
name = "2022_01_28_biasmap_noa" 
name = "2022_02_01_biasmap_noa"
name = "2022_01_27_powermap_noa"
xval, yval, cval, zval, data, rf, power = loader.load_hypercube(path, name)
#data[x, y, sd, gate] <- biasmaps
# data = data*1e9 + .2 - .005
# data = data * 1e9 +21.2
data = data * 1e9 + 1.2
# data = rf

save = True
show = True
savename = name + ""
normalize = False
contours = False
conduction = False
zerooffset = False
zmin = -40
zmax = -20
zmin = -10
zmax = 10

smooth = True
smoothval = 2

mapcolor = 'plasma'#'seismic'
mapcolor = 'seismic'
# mapcolor = 'Blues_r'


if normalize:
        savename = savename + "_norm"
if contours:
        savename = savename + "_rfcont"

# print(cval.shape)
# print(zval.shape)
clen = cval[:].size
zlen = zval[:].size

print("manipulating")

if zerooffset:
        offset = np.mean(data[:,:,0,:])
        data = data - offset

if smooth:
        for i in range(clen):
                for p in range(zlen):
                        data[:,:,i, p] = applyGauss(data[:,:,i, p], smoothval)

if conduction:
        derstep = np.abs(cval[0] - cval[1])
        data = np.gradient(data, derstep, axis = 2)
        savename = savename + "conduction"

xlow = 0
xhigh = 0
ylow = 0
yhigh = 0

print("building")

axs, figs = makemapgrid(zlen, clen, scale = 1)
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

if normalize:
    zmax = np.max(data)
    zmin = np.min(data)


print("populating")

for i in range(zlen):
        for p in range(clen):

                d = data[:,:,p,zlen - 1 - i]
                rfd = rf[:,:,p,i]
                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)

                ax[i][p].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[i][p].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
                tbias = "Gate = " + str(np.around( zval[i], 3) ) + '\nSD = ' + str(np.around( cval[p], 3) )
                # ax[i][p].text(.06, .95, tbias , color = 'k', horizontalalignment = 'left', verticalalignment = 'top', transform=ax[i][p].transAxes)

if save:
        pl.savefig("C:/Jed/Figures/" + savename + ".png")
if show:
        pl.show()