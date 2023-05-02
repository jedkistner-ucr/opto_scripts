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
from curves import simplepower, exponential, twobody

name = "CPS_2021_03_18_1"
name = "CPS_2021_11_15_18_filt"
name = "CPS_2021_11_12_42"
name = "CPS_2021_12_09_37"
name = "CPS_2022_01_26_70"
name = "2022_02_03_12_alphamap"
# name = "2022_02_03_12_errormap"

name = "2022_02_03_12_powercube"
name = "CPS_2022_01_26_53" # high power cut thing
name = "CPS_2022_01_26_55" # low power cut thing
# name = "CPS_2022_02_01_6"

path = "C:/QMO/Built"

# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
# data = data[:,:,38]
info = make_dict(info)

xlabel = "Source Drain (V)"
ylabel = 'Power (prop. mW)'

# data = data * -1
data = data * 1e9       # calibrates to real values

fitfunc = twobody

# data = np.flip(data, 1)
# data = data * -1

# xval = np.flip(xval)
# data = np.flip(data, 1)

save = True
show = True

usepower = True

#Parameters
mapcolor_name = 'plasma' #Color of map
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 50
gateCuts = []
sourceCuts = []
plotgline = []

#Custom range for source drain slices
fitder = "y"
gatemin = -2
gatemax = 8
sdmin = .7
sdmax = 1.5
blackgatemin = 0
blackgatemax = 0
hlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 300
smooth = True
smoothval = 0
postsmooth = True
smoothpostval = 0
linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = "y"
offset = 0.17
zeropoweroffset = True
setzero = True

#trims data if nonzero
xlimits = [0]
ylimits = [0]
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

# if usepower:
#     for i in range(ylen):
#         yval[i] = np.mean(power[i,:])

for i in range(xlen):
    for p in range(ylen):
        if data[p,i] == np.inf:
            data[p,i] = 0

if poweraxis == "x":
    pval = []
    for i in range(xlen):
        pval.append(np.mean(power[:,i]))
    pval = np.asfarray(pval)
    xval = pval
    if zeropoweroffset:
        offset = np.mean(data[:,0])
        data = data - offset
if poweraxis == "y":
    pval = []
    for i in range(ylen):
        pval.append(np.mean(power[i,:]))
    pval = np.asfarray(pval)
    yval = pval
    yval = yval - yval[0]

    if zeropoweroffset:
        offset1 = np.mean(data[0,:])
        data = data - offset1 - offset
    if setzero:
        data[0,:] = np.zeros(data[0,:].shape)

data = data * -1

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
    data = applyGauss(data, smoothval)

# if smooth:
#     for i in range(ylen):
#         data[i,:] = applyGauss(data[i,:], smoothval)

if linearsmooth:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmoothval)

data0 = data

if postsmooth:
    data = applyGauss(data, smoothpostval)

if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)

xlen = xval.size
ylen = yval.size

# Makes data absolute
if xabs:
    xdata = np.abs(xdata)
if yabs:
    ydata = np.abs(ydata)

if der == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    data = data * -1
    # data, xdata = np.gradient(data)
elif der == 'x_log':
    data = np.log(np.abs(data))
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data)
    # data = data * -1
elif der == 'xx':
    ydelta, data = np.gradient(data)
    ydelta, data = np.gradient(data)
elif der == 'y': 
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)
elif der == 'yy': 
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)
    data, xdelta = np.gradient(data, derstep)
elif der == "xy":
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data1 = np.gradient(data, derstep)
    data1 = data1 * -1
    derstep = np.abs(yval[0] - yval[1])
    data2, xdelta = np.gradient(data, derstep)
    # data2 = data2 * -1
elif der == "spec_1":
    data1 = data - 3.35
    ydelta, data = np.gradient(data)
    data = -data/data1
elif der == "spec_2":
    data1 = data - 3.35
    ydelta, data = np.gradient(data)
    data = -data/data1
    ydelta, data = np.gradient(data)
elif der == "spec_3":
    # diodeval = [0.03906786, 1.03426762, 0.10662857]
    diodeval = [0.03543213, 1.07187921, 0.13021503]
    # diodeval = [.03035739, 1.119459, .1256525]
    diodeval = [0.04680082, 0.98823503, 0.1204424 ]  #1/8/2021
    lineval = [0.01718924, 0.19842941]
    crossover = 4
    scalemod = np.zeros(yval.shape)
    scalemod = diode(yval, diodeval[0],diodeval[1],diodeval[2])
    linscalemod = np.zeros(yval.shape)
    linscalemod = line(yval, lineval[0], lineval[1])
    for i in range(ylen):
        # if yval[i] < crossover:
        data[i, :] = (1 / (scalemod[ylen - i - 1])) * data[i, :]
        # else:
            # data[i, :] = (1 / (linscalemod[ylen - i - 1])) * data[i, :]
elif der == "spec_3derx":
    diodeval = [0.03906786, 1.03426762, 0.10662857]
    # diodeval = [.03035739, 1.119459, .1256525]
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data)
    data = data * -1
    scalemod = np.zeros(yval.shape)
    scalemod = diode(yval, diodeval[0],diodeval[1],diodeval[2])
    for i in range(ylen):
        data[i, :] = (1 / (scalemod[ylen - i - 1])) * data[i, :]
elif der == "spec_3log":
    diodeval = [0.03906786, 1.03426762, 0.10662857]
    # diodeval = [.03035739, 1.119459, .1256525]
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data)
    data = data * -1
    scalemod = np.zeros(yval.shape)
    scalemod = diode(yval, diodeval[0],diodeval[1],diodeval[2])
    for i in range(ylen):
        data[i, :] = (1 / (scalemod[ylen - i - 1])) * data[i, :]
    # data = np.log((data))
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data)
elif der == "ideal":
    data = np.abs(data)
    data = np.log(data)
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    q = 1.602e-19 #fundamental electron charge in coloumbs
    temperature = 300 #temperature in Kelvin
    kt = 8.617333262e-5 * temperature  #KT in electron volts
    data = 1/data

elif der == "spec_34":
    diodeval = [0.04680082, 0.98823503, 0.1204424 ]  #1/8/2021
    scalemod = np.zeros(yval.shape)
    scalemod = diode(yval, diodeval[0],diodeval[1],diodeval[2])
    linscalemod = np.zeros(yval.shape)
    for i in range(ylen):
        data[i, :] = (1 / (scalemod[ylen - i - 1])) * data[i, :]
    data = np.abs(data)
    data = np.log(data)
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    data = data * 0.0256875
    data = np.reciprocal(data)
    
if xlog:
    data = np.abs(data)

mapdata = np.zeros(data.shape)
mapdata[:,:] = data[:,:]

############################################# DATA END STUFF

#Loads and builds axes
axs, axcbs, axscbs, figs = makeaxes(6, cbar = True, scalebar=True, scale = 1)
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

if gateCuts != []:
    clrgate = color_meline(colormap_name, len(gateCuts)  )
if sourceCuts != []:
    clrsd = color_meline(colormap_name, len(sourceCuts)  )


gatenorm = mpl.colors.Normalize(vmin = np.max(gateCuts), vmax = np.min(gateCuts))
sdnorm = mpl.colors.Normalize(vmin = np.max(sourceCuts), vmax = np.min(sourceCuts))
mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
mpl.colorbar.ColorbarBase(axsc[1], cmap = pl.get_cmap(colormap_name), norm = sdnorm, orientation='horizontal')

ymin1 = 10000
ymax1 = -10000
ymin2 = 10000
ymax2 = -10000

if der == "xy":
    data = data1

gateindex = [np.searchsorted(yval, gatemin), np.searchsorted(yval, gatemax)]
sdindex = [np.searchsorted(xval, sdmin), np.searchsorted(xval, sdmax)]

clr = color_meline(colormap_name, sdindex[1] - sdindex[0])

# data = data + 1.61
# data = np.abs(data)

vals = []

if fitder == "y":
    x = yval[gateindex[0]:gateindex[1]]
    cc = 0
    for i in range(sdindex[0], sdindex[1], 1):
        y = data[gateindex[0]:gateindex[1], i]

        ax[1].plot(x, y, c = clr[cc], linewidth = 1)
        # ax[2].plot(x, y, c = clr[i], linewidth = 1, linestyle = ":")

        try:
            # par, pcov = curve_fit(simplepower, x, y, maxfev = 3200)
            # ax[2].plot(x, simplepower(x, *par), c = clr[i], alpha = 0.5)
            par, pcov = curve_fit(fitfunc, x, y, maxfev = 3200)
            error = np.sum( (np.abs(y) - np.abs(fitfunc(x, *par))) )
            e0 = np.abs(pcov[0,0])
            # ax[5].scatter(xval[i], error, c = clr[cc], s = 5)
            ax[5].scatter(xval[i], e0, c = clr[cc], s = 5)
            # if error < 1.37:
            if e0 < 10:
                ax[2].plot(x, fitfunc(x, *par), c = clr[cc], alpha = 0.5)
                ax[3].scatter(xval[i], par[0], c = clr[cc], s = 5)
                ax[4].scatter(xval[i], par[1], c = clr[cc], s = 5)

            else:
                par = [0,0]
                # ax[2].plot(x, fitfunc(x, *par), c = clr[i], alpha = 0.5)
                ax[3].scatter(xval[i], par[0], c = clr[cc], s = 5)
                ax[4].scatter(xval[i], par[1], c = clr[cc], s = 5)
            vals.append("Vsd = " + str(xval[i]) + "  Voffset = " + str(par[0]))
            
        except RuntimeError:
            None

        ax[1].plot(x, y, c = clr[cc], linewidth = 1)
        cc += 1
        



if gatemin > np.min(yval):
    ax[0].axhline(gatemin, linestyle = ":", c = 'k', linewidth = 2)
if gatemax < np.max(yval):
    ax[0].axhline(gatemax, linestyle = ":", c = 'k', linewidth = 2)
    

if sdmin > np.min(xval):
    ax[0].axvline(sdmin, linestyle = ":", c = 'k', linewidth = 2)
if sdmax < np.max(xval):
    ax[0].axvline(sdmax, linestyle = ":", c = 'k', linewidth = 2)
    

ymaxrange1 = ymax1 - ymin1
ymaxrange2 = ymax2 - ymin2
margin1 = np.abs(ymaxrange1) * .04
margin2 = np.abs(ymaxrange2) * .05


for i in hlines:
    ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = 1)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits == [0]:
    if fitder == "y":
        ax[1].set_xlim( yval[gateindex[0]], yval[gateindex[1]-1] )
        ax[2].set_xlim( yval[gateindex[0]], yval[gateindex[1]-1] )
    else:
        ax[1].set_xlim(xval[sdindex[0]], xval[sdindex[1]-1])
        ax[2].set_xlim(xval[sdindex[0]], xval[sdindex[1]-1])
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
savename = name + "_fit_" + str(today)

ax[0].set_title("Map")
ax[1].set_title("Raw Data")
ax[2].set_title("Fits")
ax[3].set_title("Par[0]")
ax[4].set_title("Par[1]")


ax[0].set_xlabel(xlabel)
ax[1].set_xlabel(xlabel)
ax[2].set_xlabel(xlabel)
ax[3].set_xlabel(xlabel)
ax[4].set_xlabel(xlabel)

for v in vals:
    print(v)

if save:
    pl.savefig("C:/QMO/generated plots/" + savename + ".png")

if show:
    pl.show()