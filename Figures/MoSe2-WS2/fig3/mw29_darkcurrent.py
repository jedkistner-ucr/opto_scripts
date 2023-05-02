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
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader



font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

mpl.rcParams['font.family'] = 'Helvetica World'

name = "CPS_2021_03_18_1" #mw29 dark

path = "C:/QMO/Built"
savename = "MW29_darkcurrent"


# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)

data = data * -1e9       # calibrates to real values
# data = data + 2.137

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\it{I}$ (nA)"
ylabel = "$\itV_{g}$ (V)"
fontsize = 22
ticksize = 22
insetfontsize = 15
linewidth = 6

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'magma' #Color of map
mapcolor_name = "viridis"
mapcolor_name = "Blues"
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
zlimit = [-2, 8]
zticks = [-2,3,8]
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

#Loads and builds axes
axs, figs = customPlot(4, 4, insetbar = True)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .8)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'w', linewidth = 2, alpha = .5)

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

ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = upper)

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[1].imshow(colordata, cmap = cmap, norm = cnorm, extent = colorxt, aspect = 'auto', origin = upper)
ax[1].set_xticks([])
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position('right')
ax[1].set_yticks(zticks, color = 'w')

if contours:
    ax[0].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[0].set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[0].set_ylim(ylimits[0], ylimits[1])


ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax[1].tick_params(axis='y', labelsize = insetfontsize, direction = 'out', color = 'k', labelcolor = 'k')
# ax[1].yaxis.label.set_color('white')
# ax[1].spines['top'].set_color('white')
# ax[1].spines['bottom'].set_color('white')
# ax[1].spines['left'].set_color('white')
# ax[1].spines['right'].set_color('white')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 5)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 5)
ax[1].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = 'k')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()