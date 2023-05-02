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

def line(x, m, b):
    return m*x+b

name = "2021_03_25_51_powercube"

path = "C:/QMO/Built"
savename = "MW29_sd_alpha_pointcheck_smoothcheck"


xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
data = data * 1e9       # calibrates to real values
dataoffset = 0
# data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
# data = data[:,:,33]

xlabel = "$V_{SD}$ (V)"
ylabel = "A  (Arb)"
zlabel = "Laser Power  (mW)"
fontsize = 17
ticksize = 13
linewidth = 2
markersize = 30

save = True
show = True

mapcolor_name = 'plasma' #Color of map
colormap_name = 'viridis_r' #Color of slices

gate = -6.5
sd = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4]
powerslice = [0, .01, .02, .03] # Go check cval after power arranging
speccolor = []

hlines = [0]
vlines = [0]

smooth = [0, .2, .5, .6, .7, .8]
linearsmooth = [0, .5, 1, 1.5, 2, 2.5, 3]
savgolsmooth = [4, 5, 6, 7, 9, 11, 13]
savgolpoly = 3

poweraxis = "z"
zeropoweroffset = True
zerofit = True
powercutoff = .14

#trims data if nonzero
xlimits = []
ylimits = []


############################################# DATA DO STUFF

xlen = xval.size
ylen = yval.size

if poweraxis == "z":
    if zeropoweroffset:
        offset_ = np.mean(data[:,:,0])
        data = data - offset_ + dataoffset
        cval = cval - np.min(cval)

# if smooth > 0:
#     for i in range(cval.size):
#         data[:,:,i] = ndimage.gaussian_filter(data[:,:,i], smooth)

powerslice = cval[:8:1]

############################################# DATA END STUFF

#Loads and builds axes
axs, figs = make_grid(len(smooth), len(sd))
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))


############################################# PLOT STUFF


gateindex = np.searchsorted(yval, gate)
powerindex = np.searchsorted(cval, powercutoff)

dmax = -np.inf
dmin = np.inf

# for i in range(len(sd)):
#     for q in range(len(linearsmooth)):

#         sdindex = np.searchsorted(xval, sd[i])

#         y = data[gateindex, sdindex, :]

#         # if savgolsmooth[q] > 1:
#         #     y = savgol_filter(y, savgolsmooth[q], savgolpoly)
#         if linearsmooth[q] > 0:
#             y = ndimage.gaussian_filter1d(y,linearsmooth[q])

#         if zerofit:
#             y = y - y[0]
#             x = x - x[0]

#         if np.max(y) > dmax:
#             dmax = np.max(y)
#         if np.min(y) < dmin:
#             dmin = np.min(y)

#         ax[i + (q*len(sd))].scatter(x[:powerindex], y[:powerindex], color = 'k')
#         # ax[i].set_title("VG,VSD,smooth = %.2f , %.2f , %.1f" % (gate, sd[i], smooth[i]), fontsize = 10)

#         try:
#             par, pcov = op.curve_fit(twobody, x[:powerindex], y[:powerindex])
#             ax[i + (q*len(sd))].plot(x[:powerindex], twobody(x[:powerindex], *par), c = 'r')
#         except:
#             print('failed')


for q in range(len(smooth)):

    data_ = np.zeros(data.shape)
    for p in range(cval.size):
        data_[:,:,p] = ndimage.gaussian_filter(data[:,:,p], smooth[q])

    for i in range(len(sd)):
    
        sdindex = np.searchsorted(xval, sd[i])

        y = data_[gateindex, sdindex, :]
        x = power[gateindex, sdindex, :]

        # if savgolsmooth[q] > 1:
        #     y = savgol_filter(y, savgolsmooth[q], savgolpoly)
        # if linearsmooth[q] > 0:
        #     y = ndimage.gaussian_filter1d(y,linearsmooth[q])

        if zerofit:
            y = y - y[0]
            x = x - x[0]

        if np.max(y) > dmax:
            dmax = np.max(y)
        if np.min(y) < dmin:
            dmin = np.min(y)

        ax[i + (q*len(sd))].scatter(x[:powerindex], y[:powerindex], color = 'k')
        ax[i + (q*len(sd))].set_title(str(sd[i]) + " - " + str(smooth[q]))
        # ax[i].set_title("VG,VSD,smooth = %.2f , %.2f , %.1f" % (gate, sd[i], smooth[i]), fontsize = 10)

        try:
            par, pcov = op.curve_fit(twobody, x[:powerindex], y[:powerindex])
            ax[i + (q*len(sd))].plot(x[:powerindex], twobody(x[:powerindex], *par), c = 'r')
        except:
            print('failed')


for i in hlines:
    for q in range(len(ax)):
        ax[q].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
for i in vlines:
    for q in range(len(ax)):
        ax[q].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits != []:
    ax[0].set_xlim(xlimits[0], xlimits[1])
if ylimits != []:
    ax[0].set_ylim(ylimits[0], ylimits[1])

# for q in range(len(ax)):
#     if dmin > 0:
#         ax[q].set_ylim(dmin * .9, dmax * 1.1)
#     else:
#         ax[q].set_ylim(dmin * 1.1, dmax * 1.1)
#     ax[q].set_xlabel(xlabel)
# ax[0].set_ylabel(ylabel)

# ax[0].tick_params(axis='x', labelsize = ticksize)
# ax[0].tick_params(axis='y', labelsize = ticksize)
# ax[0].set_xlabel(xlabel, fontsize = fontsize)
# ax[0].set_ylabel(ylabel, fontsize = fontsize)

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/meta/" + savename + ".png")

if show:
    pl.show()