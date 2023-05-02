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
alpha, beta, cmap, alphaerror, betaerror, cerror, noise, yrange, xval, yval, cval, rfi = loader.load_twobodyfit_cube(path, name)
#data[x, y, alpha] <- 
# alpha = yrange - noise

spot = [xval.size // 2, yval.size // 2]
spot = [-58.5, -36]
spot = [-57.75, -38.45]
# spot = [-52, -33]
displayindex = 18
selectmaps = [0, 1.05, 1.2, 1.6]
selectmaps = []
upperlimit = 1000
noisefloor = 1

save = True
show = True
savename = name + "_fit"
normalize = False
contours = False
zmin = 0
zmax = 10
# zmin = -.1
# zmax = .1

smooth = 1

xlim = [-65,xval[-1]]
ymin = xval[-1] + 65
ylim = [yval[-1] - ymin, yval[-1]]

mapcolor = 'plasma'#'seismic'
mapcolor = 'Blues'
mapcolor = 'Reds_r'
mapcolor = 'seismic'

clen = cval.size
xlen = xval.size
ylen = yval.size

print("manipulating")

if upperlimit > 0:
        for q in range(clen):
                for i in range(xlen):
                        for p in range(ylen):
                                if alpha[p,i,q] > upperlimit:
                                        alpha[p,i,q] = 0
if noisefloor > 0:
        for q in range(clen):
                for i in range(xlen):
                        for p in range(ylen):
                                if noisefloor * noise[p,i,q] > yrange[p,i,q]:
                                        alpha[p,i,q] = 0

if smooth > 0:
        for i in range(clen):
                alpha[:,:,i] = applyGauss(alpha[:,:,i], smooth)

spot[0] = np.searchsorted(xval, spot[0])
spot[1] = np.searchsorted(yval, spot[1])
if displayindex < 0:
        displayindex = clen + displayindex

print("building")

cols, cols_ = clen, clen
if len(selectmaps) > 0:
        cols_ = len(selectmaps)

axs, figs = makemapgrid(1, cols_)
fig = pl.figure(figsize=figs)
ax = []
for h in axs:
        for a in h:
                ax.append(fig.add_axes(a))

axs, figs = makeaxes(1)
fig1 = pl.figure(figsize=figs)
ax1 = fig1.add_axes(axs[0])

axs, figs = makeaxes(1)
fig2 = pl.figure(figsize=figs)
ax2 = fig2.add_axes(axs[0])

for i in range(len(ax)-1):
        ax[i+1].set_yticks([])

xt = get_extent(xval, yval, inverty= True)

if normalize:
    zmax = np.max(alpha)
    zmin = np.min(alpha)

print("populating")

allowedind = [0]*cols
if len(selectmaps) > 0:
        for i in range(len(selectmaps)):
                allowedind[np.searchsorted(cval, selectmaps[i])] = 1
else:
        allowedind = [1]*cols

pno = 0
for i in range(cols):
        if allowedind[i] == 1:

                d = alpha[:,:,i]
                # rfd = rfi[:,:,i]

                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)

                ax[pno].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[pno].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
                if i == displayindex:
                        ax1.imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                        if contours:
                                ax1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

                tbias = str(i) + " = " + str(np.around( cval[i], 3) )
                ax[pno].set_title(tbias)

                if xlim != (0):
                        ax[pno].set_xlim(xlim)
                if ylim != (0):
                        ax[pno].set_ylim(ylim)
                pno += 1
        

x = cval
y = np.zeros(alpha[spot[1], spot[0], :].shape)

for n in range(5):
        for m in range(5):
                y = y + alpha[spot[1] -2 + n, spot[0] -1 + m, :]
                ax1.scatter(xval[spot[0] -1 + m], yval[spot[1] -1 + n], color = 'k', alpha = .5)
y = y * (1/25)
ax2.plot(x, y)

if save:
        fig.savefig("C:/QMO/generated plots/" + savename + ".png")
if show:
        pl.show()