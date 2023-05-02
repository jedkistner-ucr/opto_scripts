'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib.font_manager as font_manager
import matplotlib as mpl
font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)
mpl.rcParams['font.family'] = 'Helvetica World'
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt, savgol_filter
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader


name = "CPS_2021_03_24_23" #mw29 initial measurement

path = "C:/QMO/Built"
savename = "MW29_pci_v_sd_1"


xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
data = data - np.mean(data[:,0:3])
data = data * 1e9       # calibrates to real values
# data = data + .009

xlabel = "$\itV_{sd}$ (V)"
zlabel = "$\itV_{g}$ (V)"
ylabel = "$\itI_{pc}$ (nA)"
fontsize = 22
ticksize = 22
insetfontsize = 15
linewidth = 6

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'PuBu' #Color of map
# mapcolor_name = "plasma_r"
# mapcolor_name = "cool"
# mapcolor_name = "viridis_r"
barcolor = 'PuBu_r'
# barcolor = 'viridis'
# barcolor = 'plasma'
# barcolor = 'cool_r'
colorpadding = 1
zerocolor = False
contours = False

invertyaxis = True
upper = "lower"

hlines = [0]
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
ylimits = [-2, 32.5]
yticks = np.arange(0, ylimits[1], 5)
# xlimits = [0]
# ylimits = [0]
zlimit = [0, -1, -2, -3, -4]
zticks = [0, -2, -4]
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

#Loads and builds axes
axs, figs = customPlot(4, 4, insetbar = True)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .8)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)


if clr == []:
    clr = color_meline(mapcolor_name, len(zlimit) + colorpadding)

x = xval
for i in range(len(zlimit)):
    zindex = np.searchsorted(yval, zlimit[i])
    y = data[zindex,:]
    if savgolsmooth > 1:
        y = savgol_filter(y, savgolsmooth, savgolpoly)

    ax[0].plot(x, y, linewidth = linewidth, c = clr[i+colorpadding])

colordata = np.transpose(np.asarray([np.linspace(zlimit[0], zlimit[-1], 1000), np.linspace(zlimit[0], zlimit[-1], 1000)]))
colorxt=[-1, 1, zlimit[0], zlimit[-1]]
ax[1].imshow(colordata, cmap = pl.get_cmap(barcolor), norm = mpl.colors.Normalize(vmin = zlimit[-1], vmax = zlimit[0]), extent = colorxt, aspect = 'auto', origin = upper)
ax[1].set_xticks([])
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position('right')
ax[1].set_yticks(zticks, color = 'b')

#Sets edges of image
if xlimits == []:
    None
else:
    ax[0].set_xlim(xlimits[0], xlimits[1])
    ax[0].set_xticks(xticks)
    ax[0].tick_params(axis='x', labelsize = ticksize)

if ylimits == []:
    None
else:
    ax[0].set_ylim(ylimits[0], ylimits[1])

ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax[1].tick_params(axis='y', labelsize = insetfontsize, direction = 'in', color = 'k', labelcolor = 'k')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 5)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 5)
ax[1].set_ylabel(zlabel, fontsize = insetfontsize, labelpad = 2, color = 'k')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()