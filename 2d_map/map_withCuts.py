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
# from curves import simplepower


name = "CPS_2021_03_18_1"
name = "CPS_2021_11_15_18_filt"
name = "CPS_2021_11_12_42"
name = "CPS_2021_12_09_37"
name = "CPS_2022_01_26_70"
name = "CPS_2022_01_26_55"
# name = "2022_02_03_12_powercube_alphamap"
# name = "2021_03_25_51_powercube_alphamap"
# name = "2021_03_25_51_powercube_filt_errormap0"
# name = "2022_02_03_12_powercube_filt_errormap0"
# name = "2022_02_03_12_powercube_alphamap"
# name1 = "2022_02_03_12_powercube_filt_alphamap_power"
# name = "2022_02_03_12_errormap"
name = "CPS_2021_03_24_23"
name = 'CPS_2021_03_18_1' #mw23 dark current

name = "CPS_2022_01_26_56"  #Ao02 tb bg

# name = "CPS_2022_02_01_16" #Ao02 bright current
# name = "CPS_2022_02_01_17" #Ao02 bright current
# name = "CPS_2022_01_26_56" #Ao02 bright current tg-bg

name = "CPS_2022_01_26_53"
name = "CPS_2022_01_26_55"

# name = "2022_02_03_15_powercube"
# name = "2022_01_26_62_powercube"
# name = "CPS_2022_02_01_16"  #Brigth current d02
# name = "CPS_2022_02_03_26_stitch"
# name = "CPS_2022_02_03_26"

# name = "CPS_2022_02_01_6"  #Dark current D02
# name = "CPS_2022_01_26_57"  #Dark current D02  Vtg and Vbg

# name1 = "CPS_2022_02_01_6"    #WM29 dark current
# name = "CPS_2022_01_26_56"
# name = "CPS_2021_03_24_23" #mw29 initial measurement

# name = "CPS_2023_03_28_19" #dark current ao07
# name = "CPS_2023_01_30_65" #dark current ao06

path = "C:/Jed/Built"

errorout = False
overlay = False
subtract = False
newgrid = False

xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
if errorout or overlay or subtract:
    xval1, yval1, data1, rf1, info1 = loader.load_map(path, name1, returnpower=False)
# xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
if overlay or subtract:
    data1 = (data1 * 1e9) - .444
# data = data - np.mean( data[20:80,20:80,0] )
# data = data[:,:,cval.size-1]

info = make_dict(info)

# xlabel = "V_{TG} (V)"
# ylabel = "V_{BG} (V)"
xlabel = ""
ylabel = ""
# ylabel = '$V_{BG} = V_{TG}$'
# ylabel = 'Power (mW)'

# data = data * -1
if overlay:
    data1 = data1 * -1
# data = (data * 1e9) + 1.74       # calibrates to real values
# data = (data * 1e9) + 0.02       # calibrates to real for mw29 dark current
data = data * 1e9 +2.8

save = True
show = True
lw = 1.5

#Parameters
mapcolor_name = 'magma' #Color of map
mapcolor_name = 'bone' #Color of map
mapcolor_name = 'seismic' #Color of map
# mapcolor_name = "PiYG"
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 50
gateCuts = []
sourceCuts = []
plotgline = [1,1,1, 1, 1]
plotsdline = [1,1,1, 1]

#Custom range for source drain slices
gatemin = -60
gatemax = 60
sdmin = -50
sdmax = 1.75
blackgatemin = 0
blackgatemax = 0
hlines = []
diagr = []
diagl = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 100
smooth = 1
postsmooth = True
smoothpostval = 0
linearsmooth = False
linearsmoothval = 0
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

poweraxis = "y"
offset =0# 2.05
zeropoweroffset = True

norm = False

#trims data if nonzero
zlog = False
xlimits = [0]
ylimits = [0]
ymaplimits = (0)
xmaplimits = (0)
zlimit = [-20,20] 
# zlimit = [0,0] 

overlimits = [-100, 100]
clevels = 50
smooth1 = True
smooth1val = 0

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

if yval[1] < yval[0]:
    yval = np.flip(yval)
    data = np.flip(data, 0)
if xval[1] < xval[0]:
    xval = np.flip(xval)
    data = np.flip(data, 1)

# for i in range(xlen):
#     for p in range(ylen):
#         if data[p,i] == np.inf or data[p,i] == -np.inf:
#             data[p,i] = 10000

if errorout:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(data1[p,i]) > 1:
                data[p,i] = 0

cutoff = 0
if cutoff > 0:
    for i in range(xlen):
        for p in range(ylen):
            if data[p,i] > cutoff or data[p,i] < 0:
                data[p,i] = 0

# lowfake = np.searchsorted(yval, -3.89)
# highfake = np.searchsorted(yval, -1)
# for i in range(ylen):
#     if i > lowfake and i < highfake:



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
    if zeropoweroffset:
        offset = np.mean(data[0,:])
        data = data - offset


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

if smooth > 0:
    data = applyGauss(data, smooth)

# if smooth:
#     for i in range(ylen):
#         data[i,:] = applyGauss(data[i,:], smoothval)

if linearsmooth:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmoothval)

data0 = data

if subtract or overlay:
    if smooth1:
        data1 = applyGauss(data1, smooth1val)

if postsmooth:
    data = applyGauss(data, smoothpostval)

if norm:
    data = data - np.min(data)
    data = data / np.max(data)

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
    # data = data * -1
    # data, xdata = np.gradient(data)
elif der == 'y_log':
    data = np.log(np.abs(data))
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)
    # data = data * -1
elif der == 'yy_log':
    data = np.log(np.abs(data))
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)
    data, xdelta = np.gradient(data, derstep)
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

elif der == "spec_type":
    data = np.log(np.abs(data) + .01)
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)

    
if xlog:
    data = np.abs(data)

if subtract:
    if newgrid:
        if xval[0] < xval1[0]:
            xmin = xval1[0]
        else:
            xmin = xval[0]
        if yval[0] < yval1[0]:
            ymin = yval1[0]
        else:
            ymin = yval[0]

        if xval[-1] > xval1[-1]:
            xmax = xval1[-1]
        else:
            xmax = xval[-1]
        if yval[-1] > yval1[-1]:
            ymax = yval1[-1]
        else:
            ymax = yval[-1]

        f1 = interp.interp2d(xval1, yval1, data1, kind = 'linear')
        f = interp.interp2d(xval, yval, data, kind = 'linear')
        xval_n = np.linspace(xmin, xmax, newbox)
        yval_n = np.linspace(ymin, ymax, newbox)

        data = f(xval_n, yval_n)
        data1 = f1(xval_n, yval_n)

        data = data - data1
        xval, xval1 = xval_n, xval_n
        yval, yval1 = yval_n, yval_n

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

if overlay:

    cntmap, cntnorm = color_memap('PuOr_r', data1, dmin = overlimits[0], dmax = overlimits[1])
    ax[0].contour(data1, cmap = cntmap, norm = cntnorm, levels = np.linspace(overlimits[0], overlimits[1], clevels), extent = xt, aspect = 'auto', origin = 'lower', alpha = 1)

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

for i in range(len(gateCuts)-1, -1, -1):
    if plotgline[i] == 1:
        gateindex = np.searchsorted(yval, gateCuts[i])
        x = xval
        if gateindex < data[ :, 0].size:
            y = data[gateindex, :]
            if gateCuts[i] > gatemin and gateCuts[i] < gatemax:
                if plotgline[i] == 1:
                    ax[1].plot(x, y, linewidth = lw, c = clrgate[i])

        #Gets min and maxes of data so that the script can scale properly
                if np.min(y) < ymin1:
                    ymin1 = np.min(y)
                if np.max(y) > ymax1:
                    ymax1 = np.max(y)

if der == "xy":
    data = data2

for i in range(len(sourceCuts)):
    # print(i)
    if plotsdline[i] == 1:
        sdindex = np.searchsorted(xval, sourceCuts[i])
        # if sdindex < xval.size:
        x = yval
        y = data[: , sdindex]
        if sourceCuts[i] > sdmin and sourceCuts[i] < sdmax: 
            ax[2].plot(x, y, linewidth = lw, c = clrsd[i])
        # print(i)
        #Gets min and maxes of data so that the script can scale properly
            if np.min(y) < ymin2:
                ymin2 = np.min(y)
            if np.max(y) > ymax2:
                ymax2 = np.max(y)
    
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
for a in diagr:
    ax[0].plot(xval, xval-a, c = 'k', linestyle = ":")
for a in diagl:
    ax[0].plot(xval, -xval - np.max(np.abs(xval))-a, c = 'k', linestyle = ":")

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
    ax[1].set_xlim(np.min(xval), np.max(xval))
    ax[2].set_xlim(np.min(yval), np.max(yval))
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

fig.suptitle(name + "  :  gaussianFilter = " + str(smooth) + "  :  " + interpbool, fontsize = 10)


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
        ax[1].set_xlabel(xlabel)
    if ylabel == "":
        ax[0].set_ylabel(info["Slow Axis Variable"])
        ax[2].set_xlabel(info["Slow Axis Variable"])
    else:
        ax[0].set_ylabel(ylabel)
        ax[2].set_xlabel(ylabel)
    ax[2].set_ylabel("Current (nA)")
    ax[1].set_ylabel("Current (nA)")

else:
    savename = savename + der

ax[0].set_title("")
ax[1].set_title("Horizontal Cuts")
ax[2].set_title("Vertical Cuts")

if save:
    pl.savefig("C:/Jed/Figures/" + savename +".png")

if show:
    pl.show()
