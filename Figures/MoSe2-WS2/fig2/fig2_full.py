'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib as mpl
# import matplotlib.pyplot as pl
# import matplotlib.colors as colors
# from matplotlib import cm
import matplotlib.font_manager as font_manager
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date
import sys

from plot_tools import *
from analysis_tools import *
from custom_maps import *
import loader



font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

mpl.rcParams['font.family'] = 'Helvetica World'

#Loads and builds axes
axs, figs = fig_2()
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
ax[0].set_xticks([])
ax[0].set_yticks([])

savename = "fig_2_full"

'''
.                                                                  total v data
'''

name = 'ao2_totalv' #mw29 dark

path = "C:/QMO/Built"


xval, yval, data = loader.load_simplemap(path, name)
xval = 1240/xval
data = np.flip(data, axis = 1)
xval = np.flip(xval)
data = data - np.min(data)
# data = np.log(np.abs(data)+.01)
normval = np.max(data)
data = (2.2 * data / np.max(data))
xlabel = "Photon Energy (eV)"
zlabel = "PL Intensity (a.u.)"
ylabel = "$\itV_{bg} = \itV_{tg}$ (V)"

fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 2

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'terrain' #Color of map
# mapcolor_name = "PuRd"
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = []
vlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 2

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = []
ylimits = []
xticks = [1.55, 1.60, 1.65, 1.70]
yticks = [-6, -4, -2, 0, 2, 4, 6]
# xlimits = [0]
# ylimits = [0]
zlimit = [-.01,2]
zticks = [0, 1, 2]
# zlimit = [-2,2] 

zlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA DO STUFF

xlen = xval.size
ylen = yval.size

# for i in range(xlen):
#     for p in range(ylen):
#         if data[p,i] == np.inf:
#             data[p,i] = 0


if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

if smooth > 1:
    data = applyGauss(data, smooth)

if linearsmooth:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmoothval)

if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)

xlen = xval.size
ylen = yval.size

if der == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    data = data * -1
if der == 'y':
    derstep = np.abs(yval[0] - yval[1])
    data, xdata = np.gradient(data, derstep)
    data = data * -1

if zlog:
    # data = np.sign(data) * np.log(np.abs(data + 1))
    # data = np.sign(data)
    data = np.sign(data) * np.log(np.abs(data)+1)
    # data = np.log(data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA END STUFF

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[1].axhline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .8)
for i in vlines:
    ax[1].axvline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .5)

if zerocolor:
    zmax, zmin = np.max(data), np.min(data)
    if np.abs(zmax) > np.abs(zmin):
        zmin = -zmax
    else:
        zmax = -zmin

if zlimit == []:
    zlimit = [np.min(data), np.max(data)]

cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
# cmap, cnorm = custom_color(zlimit[0], zlimit[1])
xt = get_extent(xval, yval, inverty= invertyaxis, invertx =  False)

ax[1].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-2].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[-2].set_xticks([])
# ax[1].yaxis.tick_right()
# ax[1].yaxis.set_label_position('right')
ax[-2].set_yticks(zticks, color = 'w')

if contours:
    ax[1].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[1].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[1].set_ylim(ylimits[0], ylimits[1])

if len(xticks) > 1:
    ax[1].set_xticks(xticks)
if len(yticks) > 1:
    ax[1].set_yticks(yticks)


ax[1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[1].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[-2].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = 'w', labelcolor = 'w')
# ax[1].yaxis.label.set_color('white')
ax[-2].spines['top'].set_color('white')
ax[-2].spines['bottom'].set_color('white')
ax[-2].spines['left'].set_color('white')
ax[-2].spines['right'].set_color('white')
ax[1].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[1].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[-2].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 12, color = 'w', rotation = 270)

'''
.                                                                  delta v data
'''


name = 'ao2_deltav' #mw29 dark

path = "C:/QMO/Built"


xval, yval, data = loader.load_simplemap(path, name)
xval = 1240/xval
data = np.flip(data, axis = 1)
xval = np.flip(xval)
data = data - np.min(data)
# data = np.log(np.abs(data)+.01)
data = (2.2 * data / normval)
xlabel = "Photon Energy (eV)"
zlabel = "PL Intensity (a.u.)"
ylabel = "$\itV_{bg} = \itV_{tg}$ (V)"
ylabel1 = "Electric Field (V/nm)"
fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 4

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'terrain' #Color of map
# mapcolor_name = "hot"
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = []
vlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 2

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = []
ylimits = []
xticks = [1.55, 1.60, 1.65, 1.70]
yticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8]
# xlimits = [0]
# ylimits = [0]
zlimit = [-.01,3.9]
zticks = [0, 1, 2, 3]
zlimit1 = [-0.31827728, 0.31827728]
zticks_1 = [-.3, -.2, -.1, 0, .1, .2, .3]
# zlimit = [-2,2] 

zlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

ax0 = ax[2].twinx()
ax0.set_ylim(zlimit1)
ax0.set_yticks(zticks_1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA DO STUFF

xlen = xval.size
ylen = yval.size

# for i in range(xlen):
#     for p in range(ylen):
#         if data[p,i] == np.inf:
#             data[p,i] = 0


if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

if smooth > 1:
    data = applyGauss(data, smooth)

if linearsmooth:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmoothval)

if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)

xlen = xval.size
ylen = yval.size

if der == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    data = data * -1
if der == 'y':
    derstep = np.abs(yval[0] - yval[1])
    data, xdata = np.gradient(data, derstep)
    data = data * -1

if zlog:
    # data = np.sign(data) * np.log(np.abs(data + 1))
    # data = np.sign(data)
    data = np.sign(data) * np.log(np.abs(data)+1)
    # data = np.log(data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[2].axhline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .8)
for i in vlines:
    ax[2].axvline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .5)

if zerocolor:
    zmax, zmin = np.max(data), np.min(data)
    if np.abs(zmax) > np.abs(zmin):
        zmin = -zmax
    else:
        zmax = -zmin

if zlimit == []:
    zlimit = [np.min(data), np.max(data)]

cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
# cmap, cnorm = custom_color(zlimit[0], zlimit[1])
xt = get_extent(xval, yval, inverty= invertyaxis, invertx =  False)

ax[2].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-1].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[-1].set_xticks([])
# ax[1].yaxis.tick_right()
# ax[1].yaxis.set_label_position('right')
ax[-1].set_yticks(zticks, color = 'w')

if contours:
    ax[2].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[2].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[2].set_ylim(ylimits[0], ylimits[1])

if len(xticks) > 1:
    ax[2].set_xticks(xticks)
if len(yticks) > 1:
    ax[2].set_yticks(yticks)

ax0.tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[2].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[2].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[-1].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = 'w', labelcolor = 'w')
# ax[1].yaxis.label.set_color('white')
ax[-1].spines['top'].set_color('white')
ax[-1].spines['bottom'].set_color('white')
ax[-1].spines['left'].set_color('white')
ax[-1].spines['right'].set_color('white')
ax[2].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[2].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax0.set_ylabel(ylabel1, fontsize = fontsize, labelpad = 10, rotation = 270)
ax[-1].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 15, color = 'w', rotation = 270)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()