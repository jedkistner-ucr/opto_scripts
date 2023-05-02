'''
Explores a powercube
'''
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
from os.path import join
import scipy.optimize as op
from scipy.signal import butter, filtfilt
from datetime import date

from plot_tools import *
from analysis_tools import *
import loader
from curves import line, twobody, twobody1, simplepower, twobody_line
#~~~~~~~~~~~~~~~~~~~~~~Loading

name = "2022_02_03_12_powercube"
name = "2023_03_24_354_powercube"
name = "2022_02_03_32_powercube"
path = "C:/QMO/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
data = data * 1e9

#~~~~~~~~~~~~~~~~~~~~~~Parameters
save = True
show = True

fitfunc = twobody1

mapder = ''
zmin = 0
zmax = 0

mapcolor = 'plasma'
linecolor = 'plasma'
symmetricMapcolor = True


sources = [-2, -1.04, 0, .75, 1]
sources = [-4, -2, 0, 1, 4]
black_gates = [-2.3, -4, -5.5, -7]
black_gates = [-4]
modnumber = 1
# source = [[.5, 1], [.5, 1], [.5, 1]]
powerdisplayvalue = 1
smoothmap = 0
smooth = 0

dabs = False
dnorm = False

#~~~~~~~~~~~~~~~~~~~~~~Manipulation

xlen = xval.size
ylen = yval.size
plen = pval.size

# pval = np.delete(pval, 1)
# data = np.delete(data, 1, 2)
# powdata = np.delete(powdata, 1, 2)
# plen -= 1

# yval = np.flip(yval)
# xval = np.flip(xval)
# data = np.flip(data, 1)
# data = np.flip(data, 0)

pval = []
for i in range(plen):
    pval.append(np.mean(powdata[:,:,i]))
pval = np.asanyarray(pval)

print("Power values: ")
print(pval)

plen = pval.size

powerdisplayindex = np.searchsorted(pval, powerdisplayvalue)
if powerdisplayindex == plen:
    powerdisplayindex -= 1
mapdata = np.zeros(data[:,:,0].shape)
mapdata[:,:] = data[:,:,powerdisplayindex]

if smoothmap > 0:
    mapdata = applyGauss(mapdata, smoothmap)

if mapder == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, mapdata = np.gradient(mapdata[:,:], derstep)

if smooth > 0:
    for i in range(plen):
        data[:,:,i] = applyGauss(data[:,:,i], smooth)

#~~~~~~~~~~~~~~~~~~~~~~Plot Generation

axs, axcbs, axscbs, figs = makeaxes(len(sources)+1, cbar = True, scalebar=True, scale = .7)
ax = []
axsc = []
fig = pl.figure(figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
for a in axscbs:
    axsc.append(fig.add_axes(a))
for a in axsc:
    a.set_xticks([])
axcb = fig.add_axes(axcbs[0])
axcb.set_xticks([])
axcb.yaxis.tick_right()
axcb.yaxis.set_label_position('right')

#~~~~~~~~~~~~~~~~~~~~~~Plot populate

if symmetricMapcolor:
    mapcolor = 'bwr'
    zmin = -np.max(mapdata)
    zmax = np.max(mapdata)

cmap, cnorm = color_memap(mapcolor, mapdata[:,:] , dmin = zmin, dmax = zmax)
xt = get_extent(xval, yval, inverty=True)

ax[0].imshow(mapdata[:,:] , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'upper')
mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

if yval[0] > yval[-1]:
    sdindex = [np.searchsorted(-yval, -sd) for sd in black_gates]
else:
    sdindex = [np.searchsorted(yval, sd) for sd in black_gates]
for sd in black_gates:
    ax[0].axhline(sd, linewidth = 1.5, c = 'k', linestyle = ":", alpha = .5)

clr = color_meline(linecolor, ylen)

for i in range(len(sources)):
    ax[i+1].set_xlabel("Power (mW)")
    ax[i+1].set_ylabel("I (nA)")
    ax[i+1].set_title("Vsd = %.2f" % (sources[i]))

    x = pval
    if xval[0] > xval[-1]:
        gateindex = np.searchsorted(-xval, -sources[i])
    else:
        gateindex = np.searchsorted(xval, sources[i])
    for p in range(ylen):
        y = data[p, gateindex,:]
        ax[i+1].plot(x, y, linewidth = 1.5, c = clr[[p]])
    if len(sdindex) > 0:
        for sd in sdindex:
            y = data[sd, gateindex,:]
            ax[i+1].plot(x, y, linewidth = 3.5, c = 'k')


    barnorm = mpl.colors.Normalize(vmin = np.max(yval), vmax = np.min(yval))
    mpl.colorbar.ColorbarBase(axsc[i], cmap = pl.get_cmap(linecolor), norm = barnorm, orientation='horizontal')
    

ax[0].set_xlabel("Vsd (V)")
ax[0].set_ylabel("Vg (V)")

if save:
    savename = name + "_y_vg"
    pl.savefig("C:/QMO/generated plots/" + savename + ".png")

if show:
    pl.show()