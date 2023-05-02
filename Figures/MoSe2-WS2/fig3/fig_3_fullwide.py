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
from custom_maps import custom_color



font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

mpl.rcParams['font.family'] = 'Helvetica World'

savename = "fig_3_full_wide"

# axs, figs = fig_3_wide()
axs, figs = fig_3_stack()
ax = []
fig = pl.figure(num = 0, figsize=figs)
for i in range(len(axs)):
    ax.append(fig.add_axes(axs[i]))
ax[0].set_xticks([])
ax[0].set_yticks([])

'''
                                                                dark current map
'''
plotno = 1
barno = 2
name = "CPS_2021_03_18_1" #mw29 dark
path = "C:/QMO/Built"

# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)

data = data * -1e9       # calibrates to real values
# data = data + 2.137

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\it{I}$ (nA)"
ylabel = "$\itV_{g}$ (V)"
fontsize = 14
ticksize = 14
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
mapcolor_name = "bone"
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
    ax[plotno].axhline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .8)
for i in vlines:
    ax[plotno].axvline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .5)

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

ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[barno].set_xticks([])
ax[barno].yaxis.tick_right()
ax[barno].yaxis.set_label_position('right')
ax[barno].set_yticks(zticks, color = 'w')

if contours:
    ax[plotno].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])


ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color('white')
ax[barno].spines['top'].set_color('white')
ax[barno].spines['bottom'].set_color('white')
ax[barno].spines['left'].set_color('white')
ax[barno].spines['right'].set_color('white')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)


'''
                                                                pci map
'''
plotno = 4
barno = 5
name = "2021_03_25_51_powercube"
name = "2021_03_25_51_powercube"
path = "C:/QMO/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
data = data[:,:,-14] - np.mean(data[:,:,0])
data = data * 1e9       # calibrates to real values

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\itI_{pc}$ (nA)"
ylabel = "$\itV_{g}$ (V)"
fontsize = 14
ticksize = 14
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
mapcolor_name = 'Blues_r'
linecolor = 'gist_rainbow_r'
# linecolor = 'copper_r'
# linecolor = 'Reds'
# linecolor = 'OrRd'
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = [ -5.5, -6, -6.5, -7]
hlines= np.linspace(-5.5, -7, 5)
vlines = []
colorpad = 0
hlcolor = color_meline(linecolor, len(hlines) + colorpad)
hlcolor = ["blue", 'deepskyblue',  'orange', 'red','k']


#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = 2

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [.4, 2.4]
ylimits = [-8, -4]
xticks = [.5, 1, 1.5, 2]
yticks = [-4, -5, -6, -7, -8]
# xlimits = [0]
# ylimits = [0]
zlimit = [-.1, .38]
zticks = [0, .2, .4]
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
    ax[plotno].axhline(hlines[i], linestyle = "dashed", c = hlcolor[i + colorpad], linewidth = 1, alpha = .8)
for i in vlines:
    ax[plotno].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)

if zerocolor:
    zmax, zmin = np.max(data), np.min(data)
    if np.abs(zmax) > np.abs(zmin):
        zmin = -zmax
    else:
        zmax = -zmin

if zlimit == []:
    zlimit = [np.min(data), np.max(data)]

# cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
# cmap, cnorm = custom_color(zlimit[0], zlimit[1])
cmap, cnorm = custom_color(zlimit[0], zlimit[1])
# cmap = pl.get_cmap(mapcolor_name)
# cnorm = mpl.colors.LogNorm(vmin = zlimit[0], vmax = zlimit[1])
xt = get_extent(xval, yval, inverty= invertyaxis, invertx =  False)

ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[barno].set_xticks([])
ax[barno].yaxis.tick_right()
ax[barno].yaxis.set_label_position('right')
ax[barno].set_yticks(zticks, color = 'w')


if contours:
    ax[plotno].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].set_xticks(xticks)
ax[plotno].set_yticks(yticks)
ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color('white')
ax[barno].spines['top'].set_color('white')
ax[barno].spines['bottom'].set_color('white')
ax[barno].spines['left'].set_color('white')
ax[barno].spines['right'].set_color('white')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)
# ax[plotno].yaxis.tick_right()
# ax[plotno].yaxis.set_label_position('right')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

'''
                                                                pci curves
'''
plotno = 3
# barno = 6
name = "2021_03_25_51_powercube"
name = "2021_03_25_51_powercube"
path = "C:/QMO/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
data = data[:,:,-14] - np.mean(data[:,:,0])
data = data * 1e9       # calibrates to real values
data = data + .008

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\itV_{g}$ (V)"
ylabel = "$\itI_{pc}$ (nA)"
fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 4

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
smooth = 1.5
savgolsmooth = 0
savgolpoly = 3

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [xval[0], xval[-1]]
xticks = [.5, 1, 1.5, 2 , 2.5]
ylimits = []
# yticks = np.arange(0, ylimits[1], 5)
yticks = []
# xlimits = [0]
# ylimits = [0]
zlimit = [ -1, -2, -3, -4]
zlimit = [ -5.5, -6, -6.5, -7]
zlimit = np.linspace(-5.5, -7, 5)
zticks = [0, -1, -2, -3, -4]
zticks = []
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
    ax[plotno].axhline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .8)
for i in vlines:
    ax[plotno].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)


if clr == []:
    clr = color_meline(mapcolor_name, len(zlimit) + colorpadding)

x = xval
for i in range(len(zlimit)):
    zindex = np.searchsorted(yval, zlimit[i])
    y = data[zindex,:]
    if savgolsmooth > 1:
        y = savgol_filter(y, savgolsmooth, savgolpoly)

    # ax[plotno].plot(x, y, linewidth = linewidth, c = clr[i+colorpadding])
    ax[plotno].plot(x, y, linewidth = linewidth, c = hlcolor[i+colorpad])

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
    ax[plotno].set_xlim(xlimits[0], xlimits[1])
    ax[plotno].set_xticks(xticks)
    ax[plotno].tick_params(axis='x', labelsize = ticksize)

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
# ax[-3].tick_params(axis='y', labelsize = insetfontsize, direction = 'in', color = 'k', labelcolor = 'k')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
# ax[-3].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = 'k')

'''
                                                                alpha map
'''
plotno = 6
barno = 7

path = "C:/QMO/Built"
name = "2021_03_25_51_powercube_twobody_28_0.5"
name1 = "2021_03_25_51_powercube_twobody_29_0.5"

alpha, beta, cmap, alphaerror, betaerror, cerror, noise, yrange, xval, yval, rfi, info = loader.load_twobodyfit(path, name)
alpha1, beta1, cmap1, alphaerror1, betaerror1, cerror1, noise1, yrange1, xval1, yval1, rfi1, info1 = loader.load_twobodyfit(path, name1)
alpha = (alpha + alpha1) * 0.5
noise = (noise + noise1) / 2
yrange = (yrange + yrange1) / 2

data = alpha/400
error = np.sqrt(alphaerror * alphaerror + betaerror * betaerror)

xlabel = "$\itV_{sd}$ (V)"
zlabel = r'$ \gamma$ / $(\beta + \alpha)$ (arb.)'
ylabel = "$\itV_{g}$ (V)"
ylabel = ""
fontsize = 14
ticksize = 14
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
# mapcolor_name = 'magma'
linecolor = 'gist_rainbow_r'
# linecolor = 'copper_r'
# linecolor = 'Reds'
# linecolor = 'OrRd'
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = [ -5.5, -6, -6.5, -7]
vlines = []
colorpad = 0
hlcolor = color_meline(linecolor, len(hlines) + colorpad)
hlcolor = ["blue", 'deepskyblue', 'orange', 'red']


#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 400
smooth = .7
errortrim = 0
upperlimit = 10000000000
noisetrim = 2
noisetrimsdcutoff = .9

linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""

#trims data if nonzero

xlimits = [.4, 2.4]
ylimits = [-8, -4]
xticks = [ .5, 1, 1.5, 2 ]
yticks = [-8, -4, 0, 4, 8]
yticks = ([-4, -5, -6, -7, -8])
ticklabels = (['','','','',''])
# yticks = []
# xlimits = [0]
# ylimits = [0]
zlimit = [0, 6]
zticks = [0, 3, 6]
# zlimit = [-2,2] 

zlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA DO STUFF

xlen = xval.size
ylen = yval.size

if noisetrim > 0:
    for i in range(xlen):
        for p in range(ylen):
            if xval[i] < noisetrimsdcutoff:
                if noisetrim * noise[p,i] > yrange[p,i]:
                    data[p,i] = 0

if errortrim > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(error[p,i]) > errortrim:
                data[p,i] = 0
            if data[p,i] < 0:
                data[p,i] = 0
            if data[p,i] == np.inf or data[p,i] == -np.inf:
                data[p,i] = 0

if upperlimit > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(data[p,i]) > upperlimit:
                data[p,i] = 0

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

# for i in range(len(hlines)):
#     ax[plotno].axhline(hlines[i], linestyle = "dashed", c = hlcolor[i + colorpad], linewidth = 1, alpha = .8)
# for i in vlines:
#     ax[plotno].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)

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

ax[plotno].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[barno].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[barno].set_xticks([])
ax[barno].yaxis.tick_right()
ax[barno].yaxis.set_label_position('right')
ax[barno].set_yticks(zticks, color = 'w')


if contours:
    ax[plotno].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[plotno].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[plotno].set_ylim(ylimits[0], ylimits[1])

ax[plotno].set_xticks(xticks)
ax[plotno].set_yticks(yticks)
ax[plotno].set_yticklabels(ticklabels)
ax[plotno].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
ax[plotno].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'w')
ax[barno].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = insetcolor, labelcolor = insetcolor)
ax[barno].yaxis.label.set_color('white')
ax[barno].spines['top'].set_color('white')
ax[barno].spines['bottom'].set_color('white')
ax[barno].spines['left'].set_color('white')
ax[barno].spines['right'].set_color('white')
ax[plotno].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[plotno].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
ax[barno].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = insetcolor)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    savename = savename
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()