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

name = 'ao2_deltav' #mw29 dark

path = "C:/QMO/Built"
savename = "a02_highfield"

xval, yval, data = loader.load_simplemap(path, name)
xval = 1240/xval
data = np.flip(data, axis = 1)
xval = np.flip(xval)
data = (data - np.min(data))+.01
# data = np.log(np.abs(data)+.01)
data = (1 * data / 7810)
xlabel = "Photon Energy (eV)"
ylabel = "PL Intensity (a.u.)"
# ylabel = "$\itV_{bg} = \itV_{tg}$ (V)"

fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 4

save = True
show = True

# mapcolor_name = 'seismic' #Color of map
mapcolor_name = 'magma' #Color of map
mapcolor_name = "hot"
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

ycut = [5, 0]
clr = ['deepskyblue', 'k']
clr = ['r', 'b']

#trims data if nonzero

xlimits = [np.min(xval), 1.614]
ylimits = [-.2, 5.6]
xticks = [1.54, 1.56, 1.58, 1.60]
yticks = [0, 1, 2, 3, 4, 5]
# xlimits = [0]
# ylimits = [0]
zlimit = [0,2]
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

#Loads and builds axes
axs, figs = fig_1()
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 3, alpha = .8)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 3, alpha = .5)



if zlimit == []:
    zlimit = [np.min(data), np.max(data)]


x = xval
for i in range(len(ycut)):
    yindex = np.searchsorted(yval, ycut[i])
    y = data[yindex, :]
    ax[0].plot(x,y, linewidth = linewidth, c = clr[i])

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

if len(xticks) > 1:
    ax[0].set_xticks(xticks)
if len(yticks) > 1:
    ax[0].set_yticks(yticks)

for i in range(len(ax)):
    ax[i].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
    ax[i].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
    ax[i].set_xlabel(xlabel, fontsize = fontsize, labelpad = 5)
    ax[i].set_ylabel(ylabel, fontsize = fontsize, labelpad = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()