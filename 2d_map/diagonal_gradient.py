'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy.signal import butter, filtfilt

from datetime import date

import sys
sys.path.append("C:/Users/Jedki/Documents/analCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map

name = "CPS_2021_11_10_10_filt"
name = "CPS_2021_11_15_18_filt"
toeV = False
poweraxis = ""

path = "C:/Users/jedki/QMOdata/Built"

# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
# data = data[:,:,20]
info = make_dict(info)

data = data * 1e9  
offset = 0
# data = np.flip(data, 1)

mapcolor_name = 'plasma' #Color of map
colormap_name = 'viridis_r' #Color of slices

save = True
show = True

xlabel = 'Bottom Gate (V)'
ylabel = 'Top Gate (V)'

gatemin = -50
gatemax = 50
blackgatemin = 0
blackgatemax = 0
slices = 50
gateCuts = []
sourceCuts = []
plotgline = []

interpolate = False
newbox = 300
smooth = False
smoothval = 0
postsmooth = True
smoothpostval = [1,1,1]
linearsmooth = False
linearsmoothval = 0     

hlines = []
diag = False

#trims data if nonzero
xlimits = [0]
ylimits = [0]
zlimit = [0.018,0.046]
zlimit1 = [-.018,0.01]
zlimit2 = [-.01,.016]

############################################# DATA DO STUFF

data = data + offset

xlen = xval.size
ylen = yval.size

if poweraxis == "x":
    pval = []
    for i in range(xlen):
        pval.append(np.mean(power[:,i]))
    pval = np.asfarray(pval)
    xval = pval
if poweraxis == "y":
    pval = []
    for i in range(ylen):
        pval.append(np.mean(power[i,:]))
    pval = np.asfarray(pval)
    yval = pval

if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

if smooth:
    data = applyGauss(data, smoothval)

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

data1 = np.zeros(data.shape)
data2 = np.zeros(data.shape)
data1[:,:] = data[:,:]
data2[:,:] = data[:,:]

#check if the data is on a square grid
xstep = np.abs(xval[0] - xval[1])
ystep = np.abs(yval[0] - yval[1])
if xstep == ystep:
    derstep = np.sqrt( ((2 * xstep)**2) + ((2 * ystep)**2) )
    dhalf = derstep / 2

    # i+1 and p+1
    for i in range(xlen):
        for p in range(ylen):
            if i == 0 and p == ylen - 1:
                data1[p,i] = data[p,i]
            elif i == xlen - 1 and p == 0:
                data1[p,i] = data[p,i]

            elif i == 0:
                rise = data[p+1, i+1] - data[p,i]
                data1[p,i] = rise/dhalf
            elif i == xlen - 1:
                rise = data[p,i] - data[p-1, i-1]
                data1[p,i] = rise/dhalf

            elif p == 0:
                rise = data[p+1, i+1] - data[p,i]
                data1[p,i] = rise/dhalf
            elif p == ylen - 1:
                rise = data[p,i] - data[p-1, i-1]
                data1[p,i] = rise/dhalf

            else:
                rise = data[p+1, i+1] - data[p-1, i-1]
                data1[p,i] = rise/derstep

    # i+1 and p-1
    for i in range(xlen):
        for p in range(ylen):
            if i == 0 and p == 0:
                data2[p,i] = data[p,i]
            elif i == xlen - 1 and p == ylen - 1:
                data2[p,i] = data[p,i]

            elif i == 0:
                rise = data[p,i] - data[p-1, i+1]
                data2[p,i] = rise/dhalf
            elif i == xlen - 1:
                rise = data[p+1, i-1] - data[p,i]
                data2[p,i] = rise/dhalf

            elif p == 0:
                rise = data[p+1, i-1] - data[p,i]
                data2[p,i] = rise/dhalf
            elif p == ylen - 1:
                rise = data[p,i] - data[p-1, i+1]
                data2[p,i] = rise/dhalf

            else:
                rise = data[p+1, i-1] - data[p-1, i+1]
                data2[p,i] = rise/derstep
            
if postsmooth:
    data = applyGauss(data, smoothpostval[0])
    data1 = applyGauss(data1, smoothpostval[1])
    data2 = applyGauss(data2, smoothpostval[2])
############################################# DATA END STUFF

#Loads and builds axes
axs, axcbs, figs = makemaps(3, scale = 1)
ax = []
axcb = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
for a in axcbs:
    axcb.append(fig.add_axes(a))

for a in axcb:
    a.set_xticks([])
    a.yaxis.tick_right()
    a.yaxis.set_label_position('right')


############################################# PLOT STUFF

xt = get_extent(xval, yval, inverty = True)

cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb[0], cmap = cmap, norm = cnorm, orientation='vertical')

cmap1, cnorm1 = color_memap(mapcolor_name, data1 , dmin = zlimit1[0], dmax = zlimit1[1])
ax[1].imshow(data1 , cmap = cmap1, norm = cnorm1, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb[1], cmap = cmap1, norm = cnorm1, orientation='vertical')

cmap2, cnorm2 = color_memap(mapcolor_name, data2 , dmin = zlimit2[0], dmax = zlimit2[1])
ax[2].imshow(data2 , cmap = cmap2, norm = cnorm2, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb[2], cmap = cmap2, norm = cnorm2, orientation='vertical')


for i in hlines:
    ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    ax[2].axhline(i, linestyle = ":", c = 'k', linewidth = 1)

if diag:
    for i in range(3):
        ax[i].plot(xval, xval, linestyle = ":", c = 'w', alpha = .5)
        ax[i].plot(xval, -xval, linestyle = ":", c = 'w', alpha = .5)
        ax[i].plot(xval, -xval - 3, linestyle = ":", c = 'w', alpha = .5)
        ax[i].plot(xval, -xval - 2.5, linestyle = ":", c = 'w', alpha = .5)

if xlimits == [0]:
    ax[0].set_xlim(np.min(xval), np.max(xval))
    ax[1].set_xlim(np.min(xval), np.max(xval))
    ax[2].set_xlim(np.min(xval), np.max(xval))
else:
    ax[0].set_xlim(xlimits[0], xlimits[1])
    ax[1].set_xlim(xlimits[0], xlimits[1])
    ax[2].set_xlim(xlimits[0], xlimits[1])


if ylimits == [0]:
    ax[0].set_ylim(np.min(yval), np.max(yval))
    ax[1].set_ylim(np.min(yval), np.max(yval))
    ax[2].set_ylim(np.min(yval), np.max(yval))
else:
    ax[0].set_ylim(ylimits[0], ylimits[1])
    ax[1].set_ylim(ylimits[0], ylimits[1])
    ax[2].set_ylim(ylimits[0], ylimits[1])


############################################# END PLOT STUFF

fig.suptitle(name , fontsize = 10)

if xlabel == "":
    ax[0].set_xlabel(info["Fast Axis Variable"])
    ax[1].set_xlabel(info["Fast Axis Variable"])
    ax[2].set_xlabel(info["Fast Axis Variable"])
else:
    ax[0].set_xlabel(xlabel)
    ax[1].set_xlabel(xlabel)
    ax[2].set_xlabel(xlabel)
    
if ylabel == "":
    ax[0].set_ylabel(info["Slow Axis Variable"])
    ax[1].set_ylabel(info["Slow Axis Variable"])
    ax[2].set_ylabel(info["Slow Axis Variable"])
else:
    ax[0].set_ylabel(ylabel)
    ax[1].set_ylabel(ylabel)
    ax[2].set_ylabel(ylabel)

ax[0].set_title("Map")
ax[1].set_title("Rising Right Gradient")
ax[2].set_title("Dropping Right Gradient")

axcb[0].set_ylabel("Photocurrent (nA)")
axcb[2].set_ylabel("Differential Photoconductance")
axcb[1].set_ylabel("Differential Transphotoconductance")

if save:
    savename = name + "_cornerGrad"
    # today = date.today()
    # savename = name + "__" + str(today)
    pl.savefig("C:/Users/jedki/QMOdata/generated plots/" + savename + ".png")

if show:
    pl.show()