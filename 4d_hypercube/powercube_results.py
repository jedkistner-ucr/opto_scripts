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

name = "2021_03_24_powermap_noa_twobody"
alpha, beta, cmap, alphaerror, betaerror, cerror, noise, yrange, xval, yval, cval, rfi = loader.load_twobodyfit(path, name)
#data[x, y, alpha] <- 

save = True
show = True
savename = name + "_fit"
normalize = True
contours = False
zmin = 0
zmax = 0

smooth = 1

mapcolor = 'plasma'#'seismic'
# mapcolor = 'seismic'


if normalize:
        savename = savename + "_norm"
if contours:
        savename = savename + "_rfcont"

# print(cval.shape)
# print(zval.shape)
clen = cval[:].size
zlen = zval[:].size

print("manipulating")


if smooth > 0:
        for i in range(clen):
                data[:,:,i] = applyGauss(data[:,:,i], smooth)

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

                d = data[:,:,p,i]
                rfd = rf[:,:,p,i]
                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)

                ax[i][p].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[i][p].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
                tbias = "Gate = " + str(np.around( zval[i], 3) ) + '\nSD = ' + str(np.around( cval[p], 3) )
                ax[i][p].text(.06, .95, tbias , color = 'k', horizontalalignment = 'left', verticalalignment = 'top', transform=ax[i][p].transAxes)

if save:
        pl.savefig("C:/QMO/generated plots/" + savename + ".png")
if show:
        pl.show()