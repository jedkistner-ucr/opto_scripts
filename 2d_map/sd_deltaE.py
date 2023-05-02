'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

# import matplotlib as mpl
# import matplotlib.pyplot as pl
# import matplotlib.colors as colors
# from matplotlib import cm
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
from builder import build_map
from curves import simplepower






name = "2022_02_04_11_sdcube"
mapindex = 0
yset = -3


xset = 0
black_sd = .75

path = "C:/QMO/Built"

xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)

xlabel = "V_{TG} (V)"
ylabel = "V_{BG} (V)"

# data = data * -1
data = data * 1e9       # calibrates to real values

save = True
show = True

#Parameters
mapcolor_name = 'plasma' #Color of map
# mapcolor_name = 'seismic' #Color of map
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 50
gateCuts = []
sourceCuts = []
plotgline = []

#Custom range for source drain slices
gatemin = -60
gatemax = 50
sdmin = -50
sdmax = 50
blackgatemin = 0
blackgatemax = 0
hlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 300
smooth = True
smoothval = 1
postsmooth = True
smoothpostval = 0
linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = ""
offset = 0
zeropoweroffset = False

#trims data if nonzero
xlimits = [0]
ylimits = [0]
ymaplimits = (0)
zlimit = [0,0] 
# zlimit = [0,3000]
# zlimit = [-.1,3.1]  

xlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

data = data + offset

############################################# DATA DO STUFF

xlen = xval.size
ylen = yval.size
clen = cval.size

if yval[1] < yval[0]:
    yval = np.flip(yval)
    data = np.flip(data, 0)
if xval[1] < xval[0]:
    xval = np.flip(xval)
    data = np.flip(data, 1)

if gateCuts == []:
    gateCuts = np.linspace(np.min(yval), np.max(yval), slices)
    plotgline = np.ones(slices)
if sourceCuts == []:
    sourceCuts = np.linspace(np.min(xval), np.max(xval), slices)
    plotsdline = np.ones(slices)

if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)


if smooth:
    for i in range(clen):
        data[:,:,i] = applyGauss(data[:,:, i], smoothval)

if linearsmooth:
    for p in range(clen):
        for i in range(ylen):
            data[i,:, p] = applyGauss(data[i,:, p], linearsmoothval)


if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)


# Makes data absolute
if xabs:
    xdata = np.abs(xdata)
if yabs:
    ydata = np.abs(ydata)

mapdata = data[:,:,mapindex]

############################################# DATA END STUFF

#Loads and builds axes
axs, axcbs, axscbs, figs = makeaxes(3, cbar = True, scalebar=True, scale = 1)
ax = []
axsc = []
fig = pl.figure(num = 0, figsize=figs)
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

############################################# PLOT STUFF

no = len(gateCuts)
if gateCuts != []:
    clr = color_meline(colormap_name, no)  

xr = xval
yr = yval
    
cmap, cnorm = color_memap(mapcolor_name, mapdata , dmin = zlimit[0], dmax = zlimit[1])
xt = get_extent(xr, yr, inverty= True)

ax[0].imshow(mapdata , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

if contours:
    cntmap, cntnorm = color_memap('Greys', mapdata)
    ax[0].contour(mapdata, cmap = cntmap, norm = cntnorm, levels = 50, extent = xt, aspect = 'auto', origin = 'lower', alpha = .3)

clrgate = color_meline(colormap_name, clen)

clrsd = color_meline(colormap_name, clen  )


gatenorm = mpl.colors.Normalize(vmin = cval[-1], vmax = cval[0])
sdnorm = mpl.colors.Normalize(vmin = cval[-1], vmax = cval[0])
mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
mpl.colorbar.ColorbarBase(axsc[1], cmap = pl.get_cmap(colormap_name), norm = sdnorm, orientation='horizontal')

ymin1 = 10000
ymax1 = -10000
ymin2 = 10000
ymax2 = -10000

black_sd_index =np.searchsorted(cval, black_sd)

yindex = np.searchsorted(yval, yset)
for i in range(clen):
    # if i > 0 and i < ylen - 1:
        xmin = xval[0] - yval[-1]
        xmax = xval[-1] - yval[0]
        x = np.zeros((yindex))
        y = np.zeros((yindex))
        for p in range(yindex):
            y[p] = data[p,yindex - p,i]
            xx = yval[yindex - p] - xval[p]
            x[p] = xx

        if i == black_sd_index:
            ax[1].plot(x, y, linewidth = 1.5, c = 'k')
        else:
            ax[1].plot(x, y, linewidth = 1.5, c = clrgate[i])

        if np.min(y) < ymin1:
            ymin1 = np.min(y)
        if np.max(y) > ymax1:
            ymax1 = np.max(y)

yval1 = np.flip(yval) - np.max(yval) + yset
ax[0].plot(xval, yval1, c = 'k')

xindex = np.searchsorted(xval, xset)
for i in range(clen):
    # if i > 0 and i < ylen - 1:
        xmin1 = xval[0] + yval[0]
        xmax1 = xval[-1] + yval[-1]
        x = np.zeros((ylen))
        y = np.zeros((ylen))
        for p in range(ylen):
            y[p] = data[p,p,i]
            x[p] = xval[p] + yval[p]

        if i == black_sd_index:
            ax[2].plot(x, y, linewidth = 1.5, c = 'k')
        else:
            ax[2].plot(x, y, linewidth = 1.5, c = clrsd[i])

        if np.min(y) < ymin2:
            ymin2 = np.min(y)
        if np.max(y) > ymax2:
            ymax2 = np.max(y)
    
# xval1 = xval + np.min(xval) - xset
# ax[0].plot(xval1, yval, c = 'k')

if gatemin > np.min(yval):
    ax[0].axhline(gatemin, linestyle = ":", c = 'k', linewidth = 2)
    ax[2].axvline(gatemin, linestyle = ":", c = 'k', linewidth = 2)
if gatemax < np.max(yval):
    ax[0].axhline(gatemax, linestyle = ":", c = 'k', linewidth = 2)
    ax[2].axvline(gatemax, linestyle = ":", c = 'k', linewidth = 2)

if sdmin > np.min(xval):
    ax[0].axvline(sdmin, linestyle = ":", c = 'k', linewidth = 2)
    ax[1].axvline(sdmin, linestyle = ":", c = 'k', linewidth = 2)
if sdmax < np.max(xval):
    ax[0].axvline(sdmax, linestyle = ":", c = 'k', linewidth = 2)
    ax[1].axvline(sdmax, linestyle = ":", c = 'k', linewidth = 2)

ymaxrange1 = ymax1 - ymin1
ymaxrange2 = ymax2 - ymin2
margin1 = np.abs(ymaxrange1) * .04
margin2 = np.abs(ymaxrange2) * .05

if xlog:
    ax[1].set_ylim( ymin1 - (ymin1 * .05) , ymax1 + (ymax1 * .05) )
else:
    ax[1].set_ylim( ymin1 - margin1 , ymax1 + margin1 )

if ylog:
    ax[2].set_ylim( ymin2 - (ymin2 * .05) , ymax2 + (ymax2 * .05) )
else:
    ax[2].set_ylim( ymin2 - margin2 , ymax2 + margin2 )

for i in hlines:
    ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = 1)

if ymaplimits == (0):
    ax[0].set_xlim(np.min(xval), np.max(xval))
    ax[0].set_ylim(np.min(yval), np.max(yval))
else:
    ax[0].set_ylim(ymaplimits)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits == [0]:
    ax[1].set_xlim(xmin, xmax)
    ax[2].set_xlim(xmin1, xmax1)
else:
    ax[1].set_xlim(xlimits[0], xlimits[1])

if ylimits != [0]:
    ax[1].set_ylim(ylimits[0], ylimits[1])
    ax[2].set_ylim(ylimits[0], ylimits[1])

#Sets log
if xlog:
    ax[1].set_yscale("log")
if ylog:
    ax[2].set_xscale("log")


if interpolate:
    interpbool = "interpolated"
else:
    interpbool = ""

fig.suptitle(name + "  :  gaussianFilter = " + str(smoothval) + "  :  " + interpbool, fontsize = 10)


today = date.today()
savename = name + "_sd_delE_" + str(today)

ax[0].set_title("Map")
ax[1].set_title("Delta E cuts")
ax[2].set_title("Total Volt Cuts")

if save:
    pl.savefig("C:/QMO/generated plots/" + savename +".png")

if show:
    pl.show()
