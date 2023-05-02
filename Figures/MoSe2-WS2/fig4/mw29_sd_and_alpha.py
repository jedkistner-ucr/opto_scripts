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
from builder import build_map
from curves import simplepower, twobody


name = "2021_03_25_51_powercube"

path = "C:/QMO/Built"
savename = "MW29_sd_alpha"

xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
data = data * 1e9       # calibrates to real values
dataoffset = 0
# data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
# data = data[:,:,33]

xlabel = "$V_{SD}$ (V)"
ylabel = "A (Arb)"
zlabel = "Laser Power (mW)"
fontsize = 22
ticksize = 22
linewidth = 4
markersize = 30

save = True
show = True

mapcolor_name = 'plasma' #Color of map
colormap_name = 'viridis_r' #Color of slices

gate = -6.5
powerslice = [0, .01, .02, .03] # Go check cval after power arranging
speccolor = []

hlines = [0]
vlines = []

fill = [0.25, 1.047, 1.447, 2.4]
fillclr = ['limegreen', 'deepskyblue', 'darkorange']

smooth = .5
postsmooth = 1
savgolsmooth = 0
savgolpoly = 3

poweraxis = "z"
zeropoweroffset = True
zerofit = True
powercutoff = .23

#trims data if nonzero
xlimits = [.25, 2.4]
ylimits = [-150, 3350]
ylimits = [-.2, 2.2]
xticks = [.5, 1, 1.5, 2]
yticks = [0, 1, 2]


############################################# DATA DO STUFF

xlen = xval.size
ylen = yval.size

if poweraxis == "z":
    if zeropoweroffset:
        offset_ = np.mean(data[:,:,0])
        data = data - offset_ + dataoffset
        cval = cval - np.min(cval)

if smooth > 0:
    for i in range(cval.size):
        data[:,:,i] = ndimage.gaussian_filter(data[:,:,i], smooth)

powerslice = cval[:8:1]

############################################# DATA END STUFF

#Loads and builds axes
axs, figs = customPlot(5, 4)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

# axs, figs = customPlot(4, 4, sidebar = True)
# ax = []
# fig = pl.figure(num = 0, figsize=figs)
# for a in axs:
#     ax.append(fig.add_axes(a))


############################################# PLOT STUFF

x = cval
gateindex = np.searchsorted(yval, gate)
powerindex = np.searchsorted(cval, powercutoff)
sd = []
alpha = []
for i in range(xlen):
    
    y = data[gateindex, i, :]

    if savgolsmooth > 1:
        y = savgol_filter(y, savgolsmooth, savgolpoly)
    if zerofit:
        y = y - y[0]
        x = x - x[0]

    try:
        par, pcov = op.curve_fit(twobody, x[:powerindex], y[:powerindex])
        sd.append(xval[i])
        alpha.append(par[0])
    except:
        print('failed')

alpha = np.asarray(alpha)
alpha = alpha / 3000
if postsmooth > 0:
    alpha = ndimage.gaussian_filter1d(alpha, postsmooth)

ax[0].plot(sd, alpha, linewidth = linewidth, c = 'k')

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 2)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits != []:
    ax[0].set_xlim(xlimits[0], xlimits[1])
if ylimits != []:
    ax[0].set_ylim(ylimits[0], ylimits[1])

ax[0].set_xticks(xticks)
ax[0].set_yticks(yticks)

for i in range(len(fill)-1):
    ax[0].axvspan(fill[i], fill[i+1]-.01, color = fillclr[i], alpha = .2)

ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()