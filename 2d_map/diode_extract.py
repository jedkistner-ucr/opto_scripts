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
from custom_maps import SI_diode

path = "C:/Jed/Built"

name = 'CPS_2021_03_18_1' #mw23 dark current
# name = "CPS_2021_03_24_23" #mw29 initial measurement
# name = '2021_03_25_51_powercube' #powercube for mw23
# name = "CPS_2022_01_26_53"
# name = "CPS_2022_01_26_55"
# name = "CPS_2022_02_01_6"  #Dark current D02

xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
# xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
# data = data[:,:,-1]

# name1 = "CPS_2022_02_01_16"  #Dark current D02
# xval1, yval1, data1, rf1, power1, info1 = loader.load_map(path, name1, returnpower=True)

xlabel = "V_{TG} (V)"
ylabel = "V_{BG} (V)"
xlabel = ""
ylabel = ""
# ylabel = '$V_{BG} = V_{TG}$'
# ylabel = 'Power (mW)'

# data = (data * 1e9) + 1.74       # calibrates to real values
data = (data * -1e9) + 0.02       # calibrates to real for mw29 dark current
# data = data * 1e9 +2.08             #calibrates to nA for mw29 inital measurement (bright)
# data = data * 1e9 + .214          #for powercube MW23
# data = data * -1e9

save = True
show = True

#fit ranges
sdrange = [0,.8]
gaterange = [-5, 5]
# sdrange = [0,3]
# gaterange = [-6.65,-2.75]
#cutoff for a good fit (i guess i should probably use error here but whatevs)
ncutoff = 1
hlines = [-1.38]

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
diagr = []
diagl = []

#simple parameters -- in the otf version these pretty much stay the way they are
interpolate = False
newbox = 100
smooth = 0
linearsmooth = 0
postSmooth = 2
der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

#trims data if nonzero
zlog = False
xlimits = [0]
ylimits = [0]
ymaplimits = (0)
xmaplimits = (0)
zlimit = [-200,200] 
zlimit = [0,300] 
zlimit = [0,16] 
# zlimit = [0,.3] 
# zlimit = [0,0]

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

if smooth > 0:
    data = applyGauss(data, smooth)

if linearsmooth > 0:
    for i in range(ylen):
        data[i,:] = applyGauss(data[i,:], linearsmooth)

data0 = data

if interpolate:
    f = interp.interp2d(xval, yval, data, kind = 'linear')
    xval = np.linspace(np.min(xval), np.max(xval), newbox)
    yval = np.linspace(np.min(yval), np.max(yval), newbox)
    data = f(xval, yval)

xlen = xval.size
ylen = yval.size

#save raw datasets for comparison
xval0 = np.zeros(xval.shape)
xval0[:] = xval[:]
yval0 = np.zeros(yval.shape)
yval0[:] = yval[:]
data0 = np.zeros(data.shape)
data0[:,:] = data[:,:]

if der == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, data = np.gradient(data, derstep)
    # data = data * -1
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
# axs, axcbs, axscbs, figs = makeaxes(4, cbar = True, scalebar=True, scale = 1)
# ax = []
# axsc = []
# fig = pl.figure(num = 0, figsize=figs)
# for a in axs:
#     ax.append(fig.add_axes(a))
# for a in axscbs:
#     axsc.append(fig.add_axes(a))
# for a in axsc:
#     a.set_xticks([])
# axcb = fig.add_axes(axscbs[0])
# axcb.set_yticks([])
# axcb.yaxis.tick_right()
# axcb.yaxis.set_label_position('right')

# axs, axcbs, axscbs, figs = makeaxes(2, cbar = True, scalebar=True, scale = 1)
# ax1 = []
# axsc1 = []
# fig1 = pl.figure(num = 1, figsize=figs)
# for a in axs:
#     ax1.append(fig1.add_axes(a))
# for a in axscbs:
#     axsc1.append(fig1.add_axes(a))
# for a in axsc1:
#     a.set_xticks([])
# axcb1 = fig1.add_axes(axcbs[0])
# axcb1.set_xticks([])
# axcb1.yaxis.tick_right()
# axcb1.yaxis.set_label_position('right')

axs, figs = SI_diode()
ax = []
axsc = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

# pl.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fits

xlabel = 'Vsd (V)'
ylabel = 'Current (nA)'

def diode(x, i0, a):
    i = i0*np.exp(x*a)
    return i

# def diode(x, i0, a, x0):
#     i = i0*np.exp((x-x0)*a)
#     return i

#get index ranges for fits
gil, gih = np.searchsorted(yval, gaterange[0]), np.searchsorted(yval, gaterange[1])
sdl, sdh = np.searchsorted(xval, sdrange[0]), np.searchsorted(xval, sdrange[1])
if sdh > xval.size:
    sdh -= 1
if gih == yval.size:
    gih -= 1

#array for fit data
nfact = np.zeros((gih-gil))
ifact = np.zeros((gih-gil))
nerror = np.zeros((gih-gil))
ierror = np.zeros((gih-gil))
rval = np.zeros((gih-gil))

clr = color_meline(colormap_name, gih-gil) 

ymin, ymax = np.inf, -np.inf

for i in range(gil, gih):
    x = xval[sdl:sdh]
    y = data[i, sdl:sdh]
    ax[0].plot(x, y, c = clr[i-gil], linewidth = 1.5)

    if np.min(y) < ymin:
        ymin = np.min(y)
    if np.max(y) > ymax:
        ymax = np.max(y)

    # x = xval[sdh:]
    # y = data[i, sdh:]
    # ax[0].plot(x, y, c = clr[i-gil], linewidth = 1.5, linestyle = ":")

    try:
        # par, pcov = curve_fit(diode, x, y, bounds = ([.034, 0, -1], [.038, np.inf, 1]), maxfev = 3200)
        par, pcov = curve_fit(diode, x, y, maxfev = 3200)
        nfact[i-gil] = par[1]
        nerror[i-gil] = np.sqrt(pcov[1,1])
        ifact[i-gil] = par[0]
        ierror[i-gil] = np.sqrt(pcov[0,0])
        
        ax[1].plot(x, diode(x, *par), c = clr[i-gil], linewidth = 1.5)

        res_squared = np.sum((y-diode(x,*par))**2)
        sum_squared = np.sum((y-np.mean(y))**2)
        r_val = 1 - (res_squared/sum_squared)
        rval[i-gil] = r_val
    except:
        rval[i-gil] = 0
        nerror[i-gil] = 0
        ierror[i-gil] = 0

#scatter plots of fit parameters
# ax[2].scatter(nfact, yval[gil:gih], color = 'k')
ax[2].errorbar(nfact, yval[gil:gih], xerr = nerror, fmt='o', color = 'k')
# ax[3].errorbar(ifact, yval[gil:gih], xerr = ierror, fmt='o', color = 'k')
ax[3].scatter(rval, yval[gil:gih], color = 'k')
ax[3].axvline(.9, linestyle = ":", c = 'b', linewidth = 2)

for h in hlines:
    ax[2].axhline(h, linestyle = ":", c = 'r', linewidth = 2)
    # ax[3].axhline(h, linestyle = ":", c = 'r', linewidth = 2)
    ax[3].axhline(h, linestyle = ":", c = 'r', linewidth = 2)
    ax[4].axhline(h, linestyle = ":", c = 'r', linewidth = 2)
    ax[5].axhline(h, linestyle = ":", c = 'r', linewidth = 2)

# gatenorm = mpl.colors.Normalize(vmin = yval[gil], vmax = yval[gih])
# mpl.colorbar.ColorbarBase(axcb, cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')
# axcb.set_xlabel('Vg (V)')

ax[0].set_title('Raw Data')
ax[1].set_title('Fits')
ax[2].set_title('Ideality')
ax[3].set_title('R value')

ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(ylabel)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(ylabel)
ax[2].set_xlabel('n factor')
ax[2].set_ylabel('Vg (V)')
# ax[3].set_xlabel('i factor')
# ax[3].set_ylabel('Vg (V)')
ax[3].set_xlabel('r')
ax[3].set_ylabel('Vg (V)')

margin = np.abs(ymax-ymin)*.1

#special ideality range
ax[2].set_xlim(0,10)

#setting ranges
ax[0].set_xlim(sdrange[0], sdrange[1])
ax[1].set_xlim(sdrange[0], sdrange[1])
ax[0].set_ylim(ymin - (margin/5), ymax + margin )
ax[1].set_ylim(ymin - (margin/5), ymax + margin )
ax[2].set_ylim(gaterange[0], gaterange[1])
ax[3].set_ylim(gaterange[0], gaterange[1])
# ax[4].set_ylim(gaterange[0], gaterange[1])

# pl.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Data Adjustment

# data = data1 * -1e9


newdata = np.zeros(data.shape)
newx = []
newy = []
newdata = []

for i in range(yval.size):
    if i >= gil and i < gih:
        if nfact[i-gil] > ncutoff:
            newxval = (xval / nfact[i-gil]) / 1
            newx.append(newxval)
            newy.append(np.full(xval.shape, yval[i]))
            newdata.append(data[i,:])
newx = np.ravel(np.asarray(newx))
newy = np.ravel(np.asarray(newy))
newdata = np.ravel(np.asarray(newdata))
nx = np.linspace(np.min(newx), np.max(newx), 200)
ny = yval[gil:gih]
X, Y = np.meshgrid(nx, ny)
# points = np.stack((newx, newy), axis = 1)
interp = LinearNDInterpolator(np.stack((newx, newy), axis = 1), newdata, fill_value = 0)
nd = interp(X, Y)

xval, yval = nx, ny
data = nd
mapdata = nd

xlen = xval.size
ylen = yval.size

cmap, cnorm = color_memap(mapcolor_name, data0 , dmin = zlimit[0], dmax = zlimit[1])
xt = get_extent(xval, yval, inverty= True)
xt0 = get_extent(xval0, yval0, inverty= True)

# m, mapdata = np.gradient(mapdata)

#new data
ax[5].imshow(mapdata , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
#old data
ax[4].imshow(data0 , cmap = cmap, norm = cnorm, extent = xt0, aspect = 'auto', origin = 'lower')
# mpl.colorbar.ColorbarBase(axcb1, cmap = cmap, norm = cnorm, orientation='vertical')






if ymaplimits == (0):
    None
else:
    ax[4].set_ylim(ymaplimits)

if xmaplimits == (0):
    None
else:
    ax[4].set_xlim(xmaplimits)

if interpolate:
    interpbool = "interpolated"
else:
    interpbool = ""

# fig.suptitle(name + "  :  gaussianFilter = " + str(smoothval) + "  :  " + interpbool, fontsize = 10)


today = date.today()
savename1 = name + "_fits"
savename2 = name + "_newdata"

ax[4].set_title("Original Data Map")
ax[5].set_title("New Data Map")

ax[4].set_xlabel('Vsd (V)')
ax[4].set_ylabel('Vg (V)')
ax[5].set_xlabel('Interlayer Voltage (V)')
ax[5].set_ylabel('Vg (V)')
# axcb1.set_ylabel('Current (nA)')


if save:
    fig.savefig("C:/Jed/Figures/" + savename1 +".png")
    # fig1.savefig("C:/Jed/Figures/" + savename2 +".png")

if show:
    pl.show()
