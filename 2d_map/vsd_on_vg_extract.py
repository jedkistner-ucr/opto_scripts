'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''


import numpy as np
from scipy import ndimage
import scipy.optimize as op
from scipy import fftpack
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import butter, filtfilt
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map
from curves import simplepower

path = "C:/QMO/Built"

name = 'CPS_2021_03_18_1' #mw23 dark current
# name = "CPS_2021_03_24_23" #mw29 initial measurement
# name = '2021_03_25_51_powercube' #powercube for mw23
# name = "CPS_2022_01_26_53"
# name = "CPS_2022_01_26_55"
# name = "CPS_2022_02_01_6"  #Dark current D02

xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
# data = data[:,:,-1]

xlabel = "V_{TG} (V)"
ylabel = "V_{BG} (V)"
xlabel = ""
ylabel = ""
# ylabel = '$V_{BG} = V_{TG}$'
# ylabel = 'Power (mW)'

# data = (data * 1e9) + 1.74       # calibrates to real values
data = (data * -1e9) + 0.02       # calibrates to real for mw29 dark current
# data = data * 1e9 +2.08             #calibrates to nA for mw29 inital measurement (bright)
# data = data * 1e9 + .212

save = True
show = True

#fit ranges
sdrange = [-5,5]
gaterange = [-3, 2]
gatemaxrange = [-1.5,.5]
gatevariancerange = [-5, -3]


#cutoff for a good fit (i guess i should probably use error here but whatevs)
ncutoff = 1

#Parameters
mapcolor_name = 'magma' #Color of map
mapcolor_name = 'bone' #Color of map
# mapcolor_name = 'seismic' #Color of map
# mapcolor_name = "plasma"
colormap_name = 'viridis_r' #Color of slices
contours = False

#How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
slices = 50
gateCuts = []
sourceCuts = []
plotgline = [1,1,1]
plotsdline = [1,1,1]

#Custom range for source drain slices
gatemin = -20
gatemax = 20
sdmin = -50
sdmax = 50
blackgatemin = 0
blackgatemax = 0
hlines = []
diagr = []
diagl = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 500
smooth = 2
linearsmooth = 0
der = "y_log"        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

#trims data if nonzero
zlog = False
xlimits = [0]
ylimits = [0]
ymaplimits = (0)
xmaplimits = (0)
zlimit = [-200,200] 
zlimit = [0,18] 

overlimits = [-100, 100]
clevels = 50
smooth1 = True
smooth1val = 0

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

cutoff = 0
if cutoff > 0:
    for i in range(xlen):
        for p in range(ylen):
            if data[p,i] > cutoff or data[p,i] < 0:
                data[p,i] = 0




if np.min(xval) == np.max(xval):
    xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
if np.min(yval) == np.max(yval):
    yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

rawdata = np.zeros(data.shape)
rawdata[:,:] = data[:,:]
xval00 = np.zeros(xval.shape)
xval00[:] = xval[:]
yval00 = np.zeros(yval.shape)
yval00[:] = yval[:]

if smooth > 0:
    data = applyGauss(data, smooth)

if linearsmooth > 0:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmooth)

if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)
    f1 = interp.interp2d(xval00, yval00, rawdata, kind = 'linear')
    xval = np.linspace(np.min(xval00), np.max(xval00), newbox)
    yval = np.linspace(np.min(yval00), np.max(yval00), newbox)
    rawdata = f1(xval, yval)

#save raw datasets for comparison
xval0 = np.zeros(xval.shape)
xval0[:] = xval[:]
yval0 = np.zeros(yval.shape)
yval0[:] = yval[:]
data0 = np.zeros(data.shape)
data0[:,:] = data[:,:]



xlen = xval.size
ylen = yval.size

if der == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    # data = data * -1
    # data, xdata = np.gradient(data)
elif der == 'y_log':
    data = np.log(np.abs(data))
    derstep = np.abs(yval[0] - yval[1])
    data, xdelta = np.gradient(data, derstep)
    rawdata = np.log(np.abs(rawdata))
    derstep = np.abs(yval[0] - yval[1])
    rawdata, xdelta = np.gradient(rawdata, derstep)
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
elif der == "ideal":
    data = np.abs(data)
    data = np.log(data)
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    q = 1.602e-19 #fundamental electron charge in coloumbs
    temperature = 300 #temperature in Kelvin
    kt = 8.617333262e-5 * temperature  #KT in electron volts
    data = 1/data

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
# for a in axscbs:
#     axsc.append(fig.add_axes(a))
# for a in axsc:
#     a.set_xticks([])
axcb = fig.add_axes(axscbs[0])
axcb.set_yticks([])
# axcb.yaxis.tick_right()
# axcb.yaxis.set_label_position('right')

axs, axcbs, axscbs, figs = makeaxes(2, cbar = True, scalebar=True, scale = 1)
ax1 = []
axsc1 = []
fig1 = pl.figure(num = 1, figsize=figs)
for a in axs:
    ax1.append(fig1.add_axes(a))
# for a in axscbs:
#     axsc1.append(fig1.add_axes(a))
# for a in axsc1:
#     a.set_xticks([])
axcb1 = fig1.add_axes(axcbs[0])
axcb1.set_xticks([])
axcb1.yaxis.tick_right()
axcb1.yaxis.set_label_position('right')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fits

xlabel = 'Vg (V)'
ylabel = 'Current (nA)'

def diode(x, i0, a):
    i = i0*np.exp(x*a)
    return i

#get index ranges for fits
gil, gih = np.searchsorted(yval, gaterange[0]), np.searchsorted(yval, gaterange[1])
gsl, gsh = np.searchsorted(yval, gatemaxrange[0]), np.searchsorted(yval, gatemaxrange[1])
gvl, gvh = np.searchsorted(yval, gatevariancerange[0]), np.searchsorted(yval, gatevariancerange[1])
sdl, sdh = np.searchsorted(xval, sdrange[0]), np.searchsorted(xval, sdrange[1])
if sdh == xval.size:
    sdh -= 1
if gih == yval.size:
    gih -= 1

clr = color_meline(colormap_name, sdh-sdl)
# sd = []
# sd_bad = []
# maxgateindex = []
# maxgatevalue = []
# maxgatevalue_bad = []
# maxgatecurrent = []
# gatevar = []
# gatemean = []

sd = np.zeros(xval.shape)
sd_bad = np.zeros(xval.shape)
sd_good = np.zeros(xval.shape)
maxgateindex = np.zeros(xval.shape)
maxgatevalue = np.zeros(xval.shape)
maxgatevalue_bad = np.zeros(xval.shape)
maxgatecurrent = np.zeros(xval.shape)
gatevar = np.zeros(xval.shape)
gatemean = np.zeros(xval.shape)

zeroshiftindex = 0
zeroshiftval = -np.inf

for i in range(sdl, sdh):
    derstep = np.abs(yval[1] - yval[0])

    y = data[:,i]
    # y = np.gradient(y, derstep)

    maxsearch = np.argmin(y[gsl:gsh]) + gsl
    gindex = np.argmax(y[gil:maxsearch])+gil

    gatemean[i]=(np.mean(rawdata[gvl:gvh, i]))
    gatevar[i]=(np.var(rawdata[gvl:gvh, i])) 
    mgatecur = y[gindex]
    sd[i]=(xval[i])

    if mgatecur < (gatemean[i]+gatevar[i]):
        maxgatevalue_bad[i]=(yval[gindex])
        sd_bad[i]=(xval[i])
    else:
        maxgatecurrent[i]=(y[gindex])
        sd_good[i]=(xval[i])
        maxgatevalue[i]=(yval[gindex])
        maxgateindex[i]=(gindex)
        if maxgatevalue[i] > zeroshiftval:
            zeroshiftval = maxgatevalue[i]
            zeroshiftindex = gindex

    x = yval[gil:gih]
    y = data[gil:gih,i]
    # y = np.gradient(y, derstep)

    ax[0].plot(x, y, c = clr[i-sdl], linewidth = 1.5)

ax[0].scatter(maxgatevalue, maxgatecurrent, color = 'k')
ax[1].scatter(sd_good, maxgatevalue, color = 'k')
ax[1].scatter(sd_bad, maxgatevalue_bad, color = 'r')
ax[1].set_ylim(np.min(maxgatevalue)*1.2, np.max(maxgatevalue)*.8)
nn = np.min(maxgatevalue)*.8
mm = np.max(maxgatevalue)*1.2
ax[2].scatter(sd, maxgatecurrent, color = 'r')
ax[2].scatter(sd, gatemean+gatevar, color = 'b')

gatenorm = mpl.colors.Normalize(vmin = xval[sdl], vmax = xval[sdh])
mpl.colorbar.ColorbarBase(axcb, cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
axcb.set_xlabel('Vsd (V)')

ax[0].set_title('Raw Data')
ax[1].set_title('Transition Voltage')
ax[2].set_title('Noise Floor and Peak Value')

ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(ylabel)

ax[1].set_xlabel('Vsd (V)')
ax[1].set_ylabel('Vg transition (V)')

ax[2].set_xlabel('Vsd (V)')
ax[2].set_ylabel('Noise (BLUE) Value (RED)')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Data Adjustment

newdata = np.zeros(data.shape)
newx = []
newy = []
newdata = np.zeros(data0.shape)

#changes maxgateindex to be the number of vozel shifts that need to occur to align the data
for i in range(maxgateindex.size):
    if maxgateindex[i] == 0:
        None
    else:
        maxgateindex[i] = maxgateindex[i] - zeroshiftindex

maxgateindex = np.abs(np.int16(maxgateindex))
for i in range(xval.size):
    if maxgateindex[i] == 0:
        newdata[:,i] = data0[:,i]
    else:
        newline = data0[:-maxgateindex[i], i]
        newdata[maxgateindex[i]:,i] = newline[:]
# pl.figure()
# pl.scatter(sd, maxgateindex)
# pl.show()

data = newdata
mapdata = newdata

xlen = xval.size
ylen = yval.size

cmap, cnorm = color_memap(mapcolor_name, data0 , dmin = zlimit[0], dmax = zlimit[1])
xt = get_extent(xval, yval, inverty= True)
xt0 = get_extent(xval0, yval0, inverty= True)

# m, mapdata = np.gradient(mapdata)

#new data
ax1[1].imshow(mapdata , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
#old data
ax1[0].imshow(data0 , cmap = cmap, norm = cnorm, extent = xt0, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb1, cmap = cmap, norm = cnorm, orientation='vertical')

if ymaplimits == (0):
    None
else:
    ax1[0].set_ylim(ymaplimits)

if xmaplimits == (0):
    None
else:
    ax1[0].set_xlim(xmaplimits)

if interpolate:
    interpbool = "interpolated"
else:
    interpbool = ""

# fig.suptitle(name + "  :  gaussianFilter = " + str(smoothval) + "  :  " + interpbool, fontsize = 10)


today = date.today()
savename1 = name + "_ntype"
savename2 = name + "_newdata_vg"

ax1[0].set_title("Original Data Map")
ax1[1].set_title("New Data Map")

ax1[0].set_xlabel('Vsd (V)')
ax1[0].set_ylabel('Vg (V)')
ax1[1].set_xlabel('Interlayer Voltage (V)')
ax1[1].set_ylabel('Vg (V)')
axcb1.set_ylabel('Current (nA)')


if save:
    fig.savefig("C:/QMO/generated plots/" + savename1 +".png")
    fig1.savefig("C:/QMO/generated plots/" + savename2 +".png")

if show:
    pl.show()
