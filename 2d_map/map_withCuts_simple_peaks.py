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
from scipy.signal import butter, filtfilt, find_peaks
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map


name = "d01_m9_5K_plstack_tg_bg"
name = "d01_m9_5K_plstack_tg_neg_bg"
name = "d02_deltav"
name = "d02_sd"
name = "d02_totalv"

poweraxis = ""

path = "C:/QMO/Built"
xval, yval, data = loader.load_simplemap(path, name)

xlabel = ""
ylabel = ''

# data1 = applyGauss(data1, 1)
# data2 = applyGauss(data2, 1)
# data3 = applyGauss(data3, 1)
# data = data1  / data3
# data = data * data2
# data = data1 / data3
# data = data3

# xval = 1240/xval

offset = 0
# data = np.flip(data, 1)
# data = data * -1

save = True
show = True

extraplot = False

peakheight = 1000

fit = False
lowval = 0
highval = 0
p0 = [2.5, 1.7, .54, .13]
bnds = ( ( 0.1 , 0.0001 , 0.1, 0 ) , ( 6 , 2 , 1 , 1) )
#Parameters
mapcolor_name = 'plasma' #Color of map
mapcolor_name = 'magma' #Color of map
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 20
gateCuts = []
sourceCuts = []
plotgline = []

speccolor_gate = []
apsfade = False
sv = 8
apf = 2
speccolor_sd = []

#Custom range for source drain slices
gatemin = -50
gatemax = 50
blackgatemin = .5
blackgatemax = .51
hlines = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 300
smooth = True
smoothval = 2
postsmooth = True
smoothpostval = 0
linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

showfreq = False
lowpass = False
order = 3
cutoff = .1 # cutoff freq / .5 * sample rate

#trims data if nonzero
xlimits = [786, 796]
ylimits = [0]
zlimit = [0,0]
xmaplimits = (0)
ymaplimits = (0)
# zlimit = [0,3000]
# zlimit = [-.1,3.1]  

zlog = False

xlog = False
ylog = False
xabs = False    # Take absolute value of the data in x or y
yabs = False

data = data + offset

############################################# DATA DO STUFF

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


if lowpass:
    b, a = butter(order, cutoff, btype='low', analog=False)
    for i in range(ylen):
        dy = data[:, i]            
        y = filtfilt(b, a, dy)
        data[:, i] = y[:]
    
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
    
mapdata = np.zeros(data.shape)
mapdata[:,:] = data[:,:]

if zlog:
    mapdata = mapdata - np.min(mapdata) + 1
    mapdata = np.log(mapdata)

############################################# DATA END STUFF

#Loads and builds axes
axs, axcbs, axscbs, figs = makeaxes(2, cbar = True, scalebar=True, scale = 1)
ax = []
axsc = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
# for a in axscbs:
    # axsc.append(fig.add_axes(a))
for a in axsc:
    a.set_xticks([])
axcb = fig.add_axes(axcbs[0])
axcb.set_xticks([])
axcb.yaxis.tick_right()
axcb.yaxis.set_label_position('right')

if extraplot:
    axs1, figs1 = customPlot(4, 3)
    fig1 = pl.figure(figsize = figs1)
    ax2 = fig1.add_axes(axs1[0])

if fit:
    axs, figs = makeaxes(2)
    ax1 = []
    fig1 = pl.figure(num = 1, figsize=figs)
    for a in axs:
        ax1.append(fig1.add_axes(a))

############################################# PLOT STUFF

no = len(gateCuts)
if gateCuts != []:
    clr = color_meline(colormap_name, no)  

xr = xval

# if ylimits == [0]:
yr = yval
# else:
#     lowindex = np.searchsorted(yval, ylimits[0])
#     highindex = np.searchsorted(yval, ylimits[1])
#     yr = yval[lowindex : highindex ]
#     mapdata = mapdata[lowindex : highindex , :]        

cmap, cnorm = color_memap(mapcolor_name, mapdata , dmin = zlimit[0], dmax = zlimit[1])
xt = get_extent(xr, yr, inverty= True, invertx= False)

ax[0].imshow(mapdata , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

if contours:
    cntmap, cntnorm = color_memap('Greys', mapdata)
    ax[0].contour(mapdata, cmap = cntmap, norm = cntnorm, levels = 50, extent = xt, aspect = 'auto', origin = 'lower', alpha = .3)


if gateCuts != []:
    if speccolor_gate == []:
        clrgate = color_meline(colormap_name, len(gateCuts)  )
    else:
        clrgate = speccolor_gate
if sourceCuts != []:
    if speccolor_sd == []:
        clrsd = color_meline(colormap_name, len(sourceCuts)  )
    else:
        clrsd = speccolor_sd


# axsc[0].set_xticks([np.min(gateCuts), np.max(gateCuts)])
# axsc[1].set_xticks([np.min(sourceCuts), np.max(sourceCuts)])
gatenorm = mpl.colors.Normalize(vmin = np.max(gateCuts), vmax = np.min(gateCuts))
sdnorm = mpl.colors.Normalize(vmin = np.max(sourceCuts), vmax = np.min(sourceCuts))
# mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
# mpl.colorbar.ColorbarBase(axsc[1], cmap = pl.get_cmap(colormap_name), norm = sdnorm, orientation='horizontal')

ymin1 = 10000
ymax1 = -10000
ymin2 = 10000
ymax2 = -10000

if der == "xy":
    data = data1

for i in range(ylen):
    gateindex = i
    x = xval
    if gateindex < data[ :, 0].size:
        y = data[gateindex, :]

        peaks = np.asarray( find_peaks(y, height = peakheight) )
        peakindex = peaks[0]
        peakval = xval[peakindex]
        yv = np.full(peakval.shape, yval[gateindex])

        ax[0].scatter(peakval, yv, s = 3, c = 'k' )
        ax[1].scatter(peakval, yv, s = 3, c = 'k' )

if der == "xy":
    data = data2


if gatemin > np.min(yval):
    ax[0].axhline(gatemin, linestyle = ":", c = 'k', linewidth = 2)
if gatemax < np.max(yval):
    ax[0].axhline(gatemax, linestyle = ":", c = 'k', linewidth = 2)

ymaxrange1 = ymax1 - ymin1
ymaxrange2 = ymax2 - ymin2
margin1 = np.abs(ymaxrange1) * .04
margin2 = np.abs(ymaxrange2) * .05

if xlog:
    ax[1].set_ylim( ymin1 - (ymin1 * .05) , ymax1 + (ymax1 * .05) )
else:
    ax[1].set_ylim( ymin1 - margin1 , ymax1 + margin1 )

for i in hlines:
    ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = 1)

#Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
if xlimits == [0]:
    ax[1].set_xlim(np.min(xval), np.max(xval))
else:
    ax[1].set_xlim(xlimits[0], xlimits[1])

if extraplot:
    ax2.set_xlim(np.min(xval), np.max(xval))

if ylimits == [0]:
    ax[1].set_ylim(np.min(yval), np.max(yval))
else:
    # ax[1].set_ylim(ylimits[0], ylimits[1])
    # if invert_y:
    #     ax[1].set_ylim(ylimits[1], ylimits[0])
    ax[1].set_ylim(ylimits[0], ylimits[1])

if ymaplimits == (0):
    None
else:
    ax[0].set_ylim(ymaplimits)
if xmaplimits == (0):
    None
else:
    ax[0].set_xlim(xmaplimits)

if xlog:
    ax[1].set_yscale("log")

if interpolate:
    interpbool = "interpolated"
else:
    interpbool = ""

fig.suptitle(name + "  :  gaussianFilter = " + str(smoothval) + "  :  " + interpbool, fontsize = 10)

# axcb.set_ylabel(zlabel)
today = date.today()
savename = name + "__" + str(today)


ax[0].set_xlabel(xlabel)
ax[1].set_xlabel(xlabel)

savename = savename + "simple"

ax[0].set_title("Map")
ax[1].set_title("Horizontal Cuts")

if save:
    pl.savefig("C:/QMO/generated plots/" + savename + ".png")

if show:
    pl.show()