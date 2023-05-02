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
from os.path import join
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

names = ["mose2_16s_p6", "ws2_16s_p6", "hj_16s_p6"]
path = "C:/QMO/Raw Data/PL/MW29"

savename = "device_pl"

pldata = []
for i in range(len(names)):
    f = join(path, names[i] + ".txt")
    pldata.append(np.loadtxt(f))

xlabel = "Photon Energy (eV)"
ylabel = "PL Intensity (a.u.)"
# ylabel = "$\itV_{bg} = \itV_{tg}$ (V)"

fontsize = 14
ticksize = 14
insetfontsize = 12
linewidth = 2
overlapwidth = 2

log = False
smooth = True
smoothval = 5
lowpass = False
order = 5
cutoff = .02
filterspike = True
percentage = 7

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

clr = ['b', 'r', 'k']
mod = [1/6525, 1/6525, 5/6525]

#trims data if nonzero

xlimits = [1.3, 2.1]
ylimits = [-.08, 1.2]
# ylimits = []
xticks = [1.3, 1.5, 1.7, 1.9, 2.1]
yticks = [0, .5, 1]
# xlimits = [0]
# ylimits = [0]

zlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA END STUFF

#Loads and builds axes
axs, figs = customPlot(3, 2)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
ax0 = ax[0].twiny()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = linewidth, alpha = .8)
for i in vlines:
    ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = linewidth, alpha = .5)


for i in range(len(pldata)):
    x = pldata[i][:,0]
    y = pldata[i][:,1]
    if filterspike and i == 0:
        y = spikefilter(y, percentage, (200,7000))
    if lowpass:
        y = filtfilt(b, a, y)
    if smooth:
        y = applyGauss(y, smoothval) * mod[i]
    else:
        y = y * mod[i]
    if log:
        y = np.abs(y)
    if i == 2:
        ax[0].plot(x, y, linewidth = overlapwidth, c = clr[i])
    else:
        ax[0].plot(x, y, linewidth = linewidth, c = clr[i])

if log:
    ax[0].set_yscale('log')
#Sets edges of image
if xlimits == []:
    None
else:
    ax[0].set_xlim(xlimits[0], xlimits[1])
    ax0.set_xlim(xlimits[0], xlimits[1])

if ylimits == []:
    None
else:
    ax[0].set_ylim(ylimits[0], ylimits[1])

if len(xticks) > 1:
    ax[0].set_xticks(xticks)
    ax[0].set_xticks(xticks, [])
    ax0.set_xticks(xticks)
if len(yticks) > 1:
    ax[0].set_yticks(yticks)


ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax0.tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
ax0.set_xlabel(xlabel, fontsize = fontsize, labelpad = 5)
ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 5)
# ax0.xaxis.set_label_position('top') 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

if save:
    pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

if show:
    pl.show()