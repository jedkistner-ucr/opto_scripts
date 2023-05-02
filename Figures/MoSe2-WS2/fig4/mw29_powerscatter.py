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
from curves import simplepower, twobody1, exponential, twobody, line

name = "2021_03_25_51_powercube"

fit = True
fitfunc = twobody
trimval = 150

maxval = True

path = "C:/QMO/Built"
savename = "MW29_powerscatter"

xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
data = data * 1e12       # calibrates to real values
cval = cval * 1000

xlabel = r'Laser Power ($\mu $W)'
ylabel = "$\itI_{PC}$ (pA)"
fontsize = 22
ticksize = 22
linewidth = 4
markersize = 60

save = True
show = True

hlines = []
vlines = []

smooth = 1
linearsmooth = 0
savgolsmooth = 0
savgolpoly = 3

gate = [-6.5]
source = [1.15]
clr = ['k']
lineclr = ['r']

poweraxis = "z"
zeropoweroffset = True
calcpower = False
zerocurve = True

xlimits = [-.015, .2]
ylimits = [-.03, .4]
xlimits = [0, trimval]
ylimits = [0, 120]
xticks = [0, 50, 100, 150]
yticks = [.1, .2, .3, .4]
yticks = [0, 50, 100]

# xticks = [0, .5, ]


############################################# DATA DO STUFF

xlen = xval.size
ylen = yval.size

if poweraxis == "z":
    if zeropoweroffset:
        offset_ = np.mean(data[:,:,0])
        data = data - offset_
        cval = cval - np.min(cval)


if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)


############################################# DATA END STUFF

#Loads and builds axes
axs, figs = customPlot(5, 4)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

if smooth > 0:
    for p in range(cval.size):
        data[:,:,p] = ndimage.gaussian_filter(data[:,:,p], smooth)

############################################# PLOT STUFF

# data = data + offset

x = cval
for i in range(len(gate)):
    gateindex = np.searchsorted(yval, gate[i])
    sourceindex = np.searchsorted(xval, source[i])

    y = data[gateindex, sourceindex, :]

    if zerocurve:
        y = y - y[0]

    if savgolsmooth > 1:
        y = savgol_filter(y, savgolsmooth, savgolpoly)

    if linearsmooth > 0:
        y= applyGauss(y, linearsmooth)

    ax[0].scatter(x, y, s = markersize, color = clr[i], alpha = 1, zorder = 2)

    if fit:
        try:
            print("Trying fit")
            ptrim = np.searchsorted(cval, trimval)
            par, pcov = curve_fit(fitfunc, x[:ptrim], y[:ptrim], maxfev = 3200)
            print("A = %f" % (par[0]))
            print("B = %f" % (par[1]))
            x_ = np.linspace(x[0], x[-1], 1000)
            ax[0].plot(x_, fitfunc(x_, *par), linewidth = linewidth, c = lineclr[i], zorder = 1)
        except RuntimeError:
            print("Fit failed")


for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 2)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 2)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits != []:
    ax[0].set_xlim(xlimits[0], xlimits[1])
if ylimits != []:
    ax[0].set_ylim(ylimits[0], ylimits[1])

ax[0].set_xticks(xticks)
ax[0].set_yticks(yticks)

# equation = r'$\itI_{PC} = \frac{B}{A} ln(1+AP)$'
# ax[0].text(0.3, 0.1, equation, fontsize = fontsize+1, transform=ax[0].transAxes)
ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)

if save:
    print("Saving")
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()
    