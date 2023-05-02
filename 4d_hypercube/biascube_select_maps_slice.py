'''
Takes a built bias hypecube and works with it
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
import matplotlib.font_manager as font_manager
from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp

from datetime import date

from plot_tools import *
from analysis_tools import *
import loader
# from builder import build_map

font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)
mpl.rcParams['font.family'] = 'Helvetica World'
fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 4

path = "C:/QMO/Built"
name = "2021_10_29_powermap_noa"  
name = "2022_01_29_powermap"  
name = "2022_01_26_biasmap_noa" 
name = "2021_03_25_biasmap_noa"
name = "2021_03_24_powermap_noa"
name = "2023_03_28_powermap_noa"
xval, yval, cval, zval, data, rf, power = loader.load_hypercube(path, name)
#data[x, y, sd, gate] <- biasmaps
# data = data*1e9 + .2 - .005
data = data * 1e9
# data = rf

index = 2
direction = "x"
spot = [xval[xval.size // 2], yval[yval.size // 2]]
# spot = [-58.5, -36]
# spot = [-57.75, -38.45]
# spot = [-57.2, -35.8]
# spot = [-53.97, -36.14]
spotsmooth = 1
spotclr = 'k'
# spot = [-56.39, -35.26]
spot2 = [-50.2, -31.6]
spot2 = [xval[xval.size // 2], yval[yval.size // 2]]
# spot2 = [-53.9, -25.96]
spotsmooth2 = 1
spotclr2 = 'b'
displayindex = 0
selectmaps = [0, 1.05, 1.2, 1.6]
selectmaps = [0.7, 1.124, 1.288, 1.549, 2.006]
selectmaps = [0.7, 1.124, 1.316, 1.549, 2.006]


save = True
show = True
savename = name +" " + str(index) + "slice"
savename = name +" " + str(index) + "extract"
normalize = False
contours = False
der = ""
zerooffset = False
zmin = -3
zmax = 3
zmin = 0
zmax = .30

smooth = 2

xlim = [-65,xval[-1]]
ymin = xval[-1] + 65
ylim = [yval[-1] - ymin, yval[-1]]
xlim = []
ylim = []

mapcolor = 'plasma'#'seismic'
mapcolor = 'Blues_r'
barcolor = 'Blues'
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

if zerooffset:
        offset = np.mean(data[:,:,0,:])
        data = data - offset

if smooth > 0:
        for i in range(clen):
                for p in range(zlen):
                        data[:,:,i, p] = applyGauss(data[:,:,i, p], smooth)

if der == "y":
        derstep = np.abs(zval[0] - zval[1])
        data = np.gradient(data, derstep, axis = 3)
        savename = savename + "_der_y"
if der == "x":
        derstep = np.abs(cval[0] - cval[1])
        data = np.gradient(data, derstep, axis = 2)
        savename = savename + "_der_x"

xlow = 0
xhigh = 0
ylow = 0
yhigh = 0

spot[0] = np.searchsorted(xval, spot[0])
spot[1] = np.searchsorted(yval, spot[1])
spot2[0] = np.searchsorted(xval, spot2[0])
spot2[1] = np.searchsorted(yval, spot2[1])
if displayindex < 0:
        if direction == "x":
                displayindex = clen + displayindex
        elif direction == "y":
                displayindex = zlen + displayindex

print("building")

if direction == "x":
        cols = clen
        cols_ = clen
if direction == "y":
        cols = zlen
        cols_ = zlen
if len(selectmaps) > 0:
        cols_ = len(selectmaps)

axs, figs = makemapgrid(1, cols_, colorbar=True)
fig = pl.figure(figsize=figs)
ax = []
for h in axs:
        for a in h:
                ax.append(fig.add_axes(a))

axs, figs = makeaxes(1)
fig1 = pl.figure(figsize=figs)
ax1 = fig1.add_axes(axs[0])

axs, figs = makeaxes(1, scale=.86)
fig2 = pl.figure(figsize=figs)
ax2 = fig2.add_axes(axs[0])

for i in range(len(ax)-1):
        ax[i+1].set_yticks([])

xt = get_extent(xval, yval, inverty= True)

if normalize:
    zmax = np.max(data)
    zmin = np.min(data)

print("populating")

# for n in range(2):
#         for m in range(2):
#                 if direction == "y":
#                         data[spot[1] -1 + n*2, spot[0] -1 + m*2, index, :] = -3
#                 elif direction == "x":
#                         data[spot[1] -1 + n*2, spot[0] -1 + m*2, :, index] = -3
allowedind = [0]*cols
if len(selectmaps) > 0:
        for i in range(len(selectmaps)):
                if direction == "x":
                        allowedind[np.searchsorted(cval, selectmaps[i])] = 1
                elif direction == "y":
                        allowedind[np.searchsorted(zval, selectmaps[i])] = 1
else:
        allowedind = [1]*cols

pno = 0
for i in range(cols):
        if allowedind[i] == 1:
                if direction == "x":
                        d = data[:,:,i,index]
                        rfd = rf[:,:,i,index]
                if direction == "y":
                        d = data[:,:,index, i]
                        rfd = rf[:,:,index, i]

                cmap, cnorm = color_memap(mapcolor, d, dmin = zmin, dmax = zmax)

                ax[pno].imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                if contours:
                        ax[pno].contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
                if i == displayindex:
                        ax1.imshow(d, cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
                        if contours:
                                ax1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

                # if direction == "x":
                #         tbias = "X = " + str(np.around( cval[i], 3) ) + '\nY = ' + str(np.around( zval[index], 3) )
                # elif direction == "y":
                #         tbias = "X = " + str(np.around( cval[index], 3) ) + '\nY = ' + str(np.around( zval[i], 3) )
                # ax[i].text(.06, .95, tbias , color = 'k', horizontalalignment = 'left', verticalalignment = 'top', transform=ax[i].transAxes)
                if direction == "x":
                        tbias = "$\itV_{sd}$ = " + str(np.around( cval[i], 3)) + "V"
                        fig.suptitle(name + "  " + str(np.around( zval[index], 3)))
                elif direction == "y":
                        tbias = 'Y = ' + str(np.around( zval[i], 3) )
                        fig.suptitle(name + "  " + str(np.around( cval[index], 3) ))
                ax[pno].set_title(tbias, fontsize = fontsize)
                if xlim != []:
                        ax[pno].set_xlim(xlim)
                if ylim != []:
                        ax[pno].set_ylim(ylim)
                ax[pno].tick_params(axis='y', labelsize = ticksize)
                ax[pno].tick_params(axis='x', labelsize = ticksize)
                ax[pno].set_xticks([-50, -60])
                ax[pno].set_xlabel("X ($\mu$m)", fontsize = fontsize)
                pno += 1
        
ax[0].set_yticks([-30, -40])
ax[0].set_ylabel("Y ($\mu$m)", fontsize = fontsize)
colordata = np.transpose(np.asarray([np.linspace(zmin, zmax, 1000), np.linspace(zmin, zmax, 1000)]))
colorxt=[-1, 1, zmin, zmax]
ax[-1].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin ='lower')
ax[-1].set_xticks([])
ax[-1].yaxis.tick_right()
ax[-1].yaxis.set_label_position('right')
ax[-1].tick_params(axis='y', labelsize = ticksize)
ax[-1].set_yticks([0, .11, .22, .33], color = 'w')
ax[-1].set_ylabel("$\it{I}$ (nA)", fontsize = fontsize)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # spot 1
if direction == "x":
        x = cval
        y = np.zeros(data[spot[1], spot[0], :,index].shape)
elif direction == "y":
        x = zval
        y = np.zeros(data[spot[1], spot[0], index, :].shape)
for n in range(5):
        for m in range(5):
                if direction == "y":
                        y = y + data[spot[1] -2 + n, spot[0] -1 + m, index, :]
                elif direction == "x":
                        y = y + data[spot[1] -2 + n, spot[0] -1 + m, :, index]
                ax1.scatter(xval[spot[0] -1 + m], yval[spot[1] -1 + n], color = spotclr, alpha = .5)
y = y * (1/25)
if spotsmooth > 0:
        y = ndimage.gaussian_filter1d(y, spotsmooth)
ax2.plot(x, y, c = spotclr, linewidth = linewidth)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # spot 2

if direction == "x":
        x = cval
        y = np.zeros(data[spot2[1], spot2[0], :,index].shape)
elif direction == "y":
        x = zval
        y = np.zeros(data[spot2[1], spot2[0], index, :].shape)
for n in range(5):
        for m in range(5):
                if direction == "y":
                        y = y + data[spot2[1] -2 + n, spot2[0] -1 + m, index, :]
                elif direction == "x":
                        y = y + data[spot2[1] -2 + n, spot2[0] -1 + m, :, index]
                ax1.scatter(xval[spot2[0] -1 + m], yval[spot2[1] -1 + n], color = spotclr2, alpha = .5)
y = y * (1/25)
if spotsmooth2 > 0:
        y = ndimage.gaussian_filter1d(y, spotsmooth2)
ax2.plot(x, y, c = spotclr2, linewidth = linewidth)

ax2.tick_params(axis='x', labelsize = ticksize)
ax2.tick_params(axis='y', labelsize = ticksize)
ax2.set_xlim(x[0], x[-1])
ax2.set_yticks([0, .2, .4, .6])
ax2.set_ylabel("$\it{I}$ (nA)", fontsize = fontsize)
ax2.set_xticks([.7, 1.4, 2.1])
ax2.set_xticks([1, 1.5, 2])
ax2.set_xlabel("$\itV_{sd}$ (V)", fontsize = fontsize)
for i in range(len(selectmaps)):
        ax2.axvline(selectmaps[i], linestyle = ":", linewidth = 2, c = 'k', alpha = .4)


if save:
        fig.savefig("C:/QMO/generated plots/" + savename + ".png")
        fig2.savefig("C:/QMO/generated plots/" + savename + ".png")
if show:
        pl.show()