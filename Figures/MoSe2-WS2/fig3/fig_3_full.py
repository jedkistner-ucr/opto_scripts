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
from datetime import date
import sys
from plot_tools import *
from analysis_tools import *
import loader



font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

mpl.rcParams['font.family'] = 'Helvetica World'

savename = "fig_3_full"

axs, figs = fig_3()
ax = []
fig = pl.figure(num = 0, figsize=figs)
for i in range(len(axs)):
    if i == 4:
        ax.append(0)
    else:
        ax.append(fig.add_axes(axs[i]))
ax[0].set_xticks([])
ax[0].set_yticks([])

'''
                                                                dark current map
'''

name = "CPS_2021_03_18_1" #mw29 dark
path = "C:/QMO/Built"

# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)

data = data * -1e9       # calibrates to real values
# data = data + 2.137

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\it{I}$ (nA)"
ylabel = "$\itV_{g}$ (V)"
fontsize = 12
ticksize = 12
insetfontsize = 10
linewidth = 4
insetcolor = 'w'

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'magma' #Color of map
mapcolor_name = "viridis"
mapcolor_name = "Blues"
mapcolor_name = "Oranges"
mapcolor_name = "Reds_r"
mapcolor_name = "Blues_r"
# mapcolor_name = "pink"
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = []
vlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 1

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [-2, 2.5]
ylimits = []
# xlimits = [0]
# ylimits = [0]
zlimit = [0, 8]
zticks = [0,4,8]
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
xt = get_extent(xval, yval, inverty= invertyaxis, invertx =  False)

ax[1].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-2].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[-2].set_xticks([])
ax[-2].yaxis.tick_right()
ax[-2].yaxis.set_label_position('right')
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


ax[1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[1].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[-2].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[-2].yaxis.label.set_color('white')
ax[-2].spines['top'].set_color('white')
ax[-2].spines['bottom'].set_color('white')
ax[-2].spines['left'].set_color('white')
ax[-2].spines['right'].set_color('white')
ax[1].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[1].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[-2].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)


'''
                                                                pci map
'''

name = "CPS_2021_03_24_23" #mw29 initial measurement
path = "C:/QMO/Built"

xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
data = data - np.mean(data[:,0:3])
data = data * 1e9       # calibrates to real values

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\itI_{pc}$ (nA)"
ylabel = "$\itV_{g}$ (V)"
fontsize = 12
ticksize = 12
insetfontsize = 10
linewidth = 4
insetcolor = 'w'

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'magma' #Color of map
# mapcolor_name = "plasma"
mapcolor_name = 'afmhot'
mapcolor_name = 'OrRd_r'
mapcolor_name = 'PuBu_r'
linecolor = 'gist_rainbow_r'
# linecolor = 'copper_r'
# linecolor = 'Reds'
# linecolor = 'OrRd'
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = [0, -1.3, -2.6, -4]
vlines = []
colorpad = 0
hlcolor = color_meline(linecolor, len(hlines) + colorpad)
hlcolor = ["blue", 'deepskyblue', 'orange', 'red']


#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 2

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [0, 2.7]
ylimits = []
xticks = [0, 1, 2]
yticks = [-8, -4, 0, 4, 8]
# xlimits = [0]
# ylimits = [0]
zlimit = [0,26]
zticks = [0, 13, 26]
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

if smooth > 0:
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

for i in range(len(hlines)):
    ax[3].axhline(hlines[i], linestyle = "dashed", c = hlcolor[i + colorpad], linewidth = 1, alpha = .8)
for i in vlines:
    ax[3].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)

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
# cmap, cnorm = custom_color(zlimit[0], zlimit[1])
# cmap = pl.get_cmap(mapcolor_name)
# cnorm = mpl.colors.LogNorm(vmin = zlimit[0], vmax = zlimit[1])
xt = get_extent(xval, yval, inverty= invertyaxis, invertx =  False)

ax[3].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[-1].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[-1].set_xticks([])
ax[-1].yaxis.tick_right()
ax[-1].yaxis.set_label_position('right')
ax[-1].set_yticks(zticks, color = 'w')


if contours:
    ax[3].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[3].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[3].set_ylim(ylimits[0], ylimits[1])

ax[3].set_xticks(xticks)
ax[3].set_yticks(yticks)
ax[3].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[3].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax[-1].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[-1].yaxis.label.set_color('white')
ax[-1].spines['top'].set_color('white')
ax[-1].spines['bottom'].set_color('white')
ax[-1].spines['left'].set_color('white')
ax[-1].spines['right'].set_color('white')
ax[3].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[3].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[-1].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)
ax[3].yaxis.tick_right()
ax[3].yaxis.set_label_position('right')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

'''
                                                                pci curves
'''

name = "CPS_2021_03_24_23" #mw29 initial measurement
path = "C:/QMO/Built"

xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
data = data - np.mean(data[:,0:3])
data = data * 1e9       # calibrates to real values
# data = data + .009

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\itV_{g}$ (V)"
ylabel = "$\itI_{pc}$ (nA)"
fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 5

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'OrRd' #Color of map
# mapcolor_name = "plasma_r"
# mapcolor_name = "cool"
# mapcolor_name = "viridis_r"
barcolor = 'OrRd_r'
# barcolor = 'viridis'
# barcolor = 'plasma'
# barcolor = 'cool_r'
colorpadding = 1
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = []
vlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 1
savgolsmooth = 0
savgolpoly = 3

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [0, 2.7]
xticks = np.arange(xlimits[0], xlimits[1], .5)
ylimits = [-.2, 32.5]
yticks = np.arange(0, ylimits[1], 5)
# xlimits = [0]
# ylimits = [0]
zlimit = [0, -1, -2, -3, -4]
zlimit = [0, -1.3, -2.6, -4]
zticks = [0, -1, -2, -3, -4]
clr = ['k', 'b', 'r']
clr = []


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

if smooth > 0:
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
    ax[2].axhline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .8)
for i in vlines:
    ax[2].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)


if clr == []:
    clr = color_meline(mapcolor_name, len(zlimit) + colorpadding)

x = xval
for i in range(len(zlimit)):
    zindex = np.searchsorted(yval, zlimit[i])
    y = data[zindex,:]
    if savgolsmooth > 1:
        y = savgol_filter(y, savgolsmooth, savgolpoly)

    # ax[2].plot(x, y, linewidth = linewidth, c = clr[i+colorpadding])
    ax[2].plot(x, y, linewidth = linewidth, c = hlcolor[i+colorpad])

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
# ax[-3].imshow(colordata, cmap = pl.get_cmap(barcolor), norm = mpl.colors.Normalize(vmin = zlimit[-1], vmax = zlimit[0]), extent = colorxt, aspect = 'auto', origin = upper)
# ax[-3].set_xticks([])
# ax[-3].yaxis.tick_right()
# ax[-3].yaxis.set_label_position('right')
# ax[-3].set_yticks(zticks, color = 'b')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[2].set_xlim(xlimits[0], xlimits[1])
    ax[2].set_xticks(xticks)
    ax[2].tick_params(axis='x', labelsize = ticksize)

if ylimits == []:
    None
else:
    ax[2].set_ylim(ylimits[0], ylimits[1])

ax[2].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[2].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
# ax[-3].tick_params(axis='y', labelsize = insetfontsize, direction = 'in', color = 'k', labelcolor = 'k')
ax[2].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[2].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
# ax[-3].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = 'k')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    savename = savename + "_" + linecolor
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()