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

from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map

name = "CPS_2022_01_26_56" #Ao02 bright current tg-bg
# name = "CPS_2022_01_26_57"  #Dark current D02  Vtg and Vbg

path = "C:/Jed/Built"



xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
info = make_dict(info)

xlabel = "V_{TG} (V)"
ylabel = "V_{BG} (V)"
xcutlabel = "V_{TG} - V_{BG} (V)"
ycutlabel = "V_{TG} + V_{BG} (V)"
# xlabel = ""
# ylabel = ""
# ylabel = '$V_{BG} = V_{TG}$'
# ylabel = 'Power (mW)'

# data = (data * 1e9) + 1.74       # calibrates to real values
# data = (data * 1e9) + 0.02       # calibrates to real for mw29 dark current
data = data * 1e9 +2.8

save = True
show = True
lw = 1

#Parameters
mapcolor_name = 'magma' #Color of map
mapcolor_name = 'bone' #Color of map
mapcolor_name = 'seismic' #Color of map
# mapcolor_name = "PiYG"
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 50

#simple parameters -- in the otf version these pretty much stay the way they are
smooth = 1
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

#trims data if nonzero
zlog = False
xlimits = [0]
ylimits = [0]
ymaplimits = (0)
xmaplimits = (0)
zlimit = [-20,20] 
# zlimit = [0,0] 

############################################# DATA DO STUFF

xlen = xval.size
ylen = yval.size

if yval[1] < yval[0]:
    yval = np.flip(yval)
    data = np.flip(data, 0)
if xval[1] < xval[0]:
    xval = np.flip(xval)
    data = np.flip(data, 1)


if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

if smooth > 0:
    data = applyGauss(data, smooth)

data0 = data

xlen = xval.size
ylen = yval.size

mapdata = np.zeros(data.shape)
mapdata[:,:] = data[:,:]

if zlog:
    mapdata = np.log(np.abs(mapdata))

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

xr = xval
yr = yval
    
cmap, cnorm = color_memap(mapcolor_name, mapdata , dmin = zlimit[0], dmax = zlimit[1])
xt = get_extent(xr, yr, inverty= True)

ax[0].imshow(mapdata , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

gatenorm = mpl.colors.Normalize(vmin = np.min(xval)+np.min(yval), vmax = np.max(xval)+np.max(yval))
sdnorm = mpl.colors.Normalize(vmin = np.min(xval)-np.max(yval), vmax = np.max(xval)-np.min(yval))
mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
mpl.colorbar.ColorbarBase(axsc[1], cmap = pl.get_cmap(colormap_name), norm = sdnorm, orientation='horizontal')

ymin1 = np.inf
ymax1 = -np.inf
ymin2 = np.inf
ymax2 = -np.inf

diaglen = xlen + ylen
halfdiaglen = diaglen // 2
clr = color_meline(colormap_name, diaglen)

centertotalv = np.mean(yval) + np.mean(xval)
totalVstep = np.abs(yval[0]-yval[1])+np.abs(xval[0]-xval[1])

for i in range(ylen):

    y = np.diagonal(data, offset = -ylen+1+i)
    xmargin = (y.size//2) * totalVstep
    x = np.linspace(centertotalv-xmargin, centertotalv+xmargin, y.size)

    ax[2].plot(x, y, linewidth = lw, c = clr[i])

    if np.min(y) < ymin1:
        ymin1 = np.min(y)
    if np.max(y) > ymax1:
        ymax1 = np.max(y)

for i in range(xlen):

    y = np.diagonal(data, offset = i)
    xmargin = (y.size//2) * totalVstep
    x = np.linspace(centertotalv-xmargin, centertotalv+xmargin, y.size)

    ax[2].plot(x, y, linewidth = lw, c = clr[i+ylen])

    if np.min(y) < ymin1:
        ymin1 = np.min(y)
    if np.max(y) > ymax1:
        ymax1 = np.max(y)

centerdeltav = np.mean(yval) - np.mean(xval)
deltaVstep = np.abs(yval[0]-yval[1])+np.abs(xval[0]-xval[1])

data = np.rot90(data)

for i in range(ylen):

    y = np.diagonal(data, offset = -ylen+1+i)
    y = np.flip(y)
    xmargin = (y.size//2) * deltaVstep
    x = np.linspace(centerdeltav-xmargin, centerdeltav+xmargin, y.size)
    
    ax[1].plot(x, y, linewidth = lw, c = clr[i])

    if np.min(y) < ymin2:
        ymin2 = np.min(y)
    if np.max(y) > ymax2:
        ymax2 = np.max(y)

for i in range(xlen):

    y = np.diagonal(data, offset = i)
    y = np.flip(y)
    xmargin = (y.size//2) * deltaVstep
    x = np.linspace(centerdeltav-xmargin, centerdeltav+xmargin, y.size)

    ax[1].plot(x, y, linewidth = lw, c = clr[i+ylen])

    if np.min(y) < ymin2:
        ymin2 = np.min(y)
    if np.max(y) > ymax2:
        ymax2 = np.max(y)

ymaxrange1 = ymax1 - ymin1
ymaxrange2 = ymax2 - ymin2
margin1 = np.abs(ymaxrange1) * .04
margin2 = np.abs(ymaxrange2) * .05


ax[1].set_ylim( ymin1 - margin1 , ymax1 + margin1 )
# ax[2].set_ylim( ymin2 - margin2 , ymax2 + margin2 )

if ymaplimits == (0):
    None
else:
    ax[0].set_ylim(ymaplimits)
if xmaplimits == (0):
    None
else:
    ax[0].set_xlim(xmaplimits)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits == [0]:
    ax[1].set_xlim(np.max(xval)-np.min(yval), np.min(xval)-np.max(yval))
    ax[2].set_xlim(np.min(xval)+np.min(yval), np.max(xval)+np.max(yval))
else:
    ax[1].set_xlim(xlimits[0], xlimits[1])

if ylimits != [0]:
    ax[1].set_ylim(ylimits[0], ylimits[1])
    ax[2].set_ylim(ylimits[0], ylimits[1])



fig.suptitle(name + "  :  gaussianFilter = " + str(smooth) + "  :  " , fontsize = 10)


today = date.today()
savename = name + "__" + str(today)

if der == "x" or der == "xy":
    savename = savename + "xder"
    ax[1].set_ylabel("Conductance (nA/V)")
    ax[2].set_ylabel("Conductance (nA/V)")

elif der == "":
    if xlabel == "":
        ax[0].set_xlabel(info["Fast Axis Variable"])
        ax[1].set_xlabel(info["Fast Axis Variable"])
    else:
        ax[0].set_xlabel(xlabel)
        ax[1].set_xlabel(xcutlabel)
    if ylabel == "":
        ax[0].set_ylabel(info["Slow Axis Variable"])
        ax[2].set_xlabel(info["Slow Axis Variable"])
    else:
        ax[0].set_ylabel(ylabel)
        ax[2].set_xlabel(ycutlabel)
    ax[2].set_ylabel("Current (nA)")
    ax[1].set_ylabel("Current (nA)")

ax[0].set_title("")
ax[1].set_title("Delta V")
ax[2].set_title("Total V")

if save:
    pl.savefig("C:/Jed/Figures/" + savename +".png")

if show:
    pl.show()
