'''
Takes a biasmap cubed with power, extracts curves at set gate and source drain
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
from curves import twobody, twobody1, line
#~~~~~~~~~~~~~~~~~~~~~~Loading

name = "2022_02_03_12_powercube"
name = "2023_03_24_352_powercube"
name = "2021_03_25_51_powercube"    #sd-gate powercube for MW23
name = '2021_03_25_51_powercube_butter_4_2.0'
path = "C:/Jed/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
data = data * 1e9
#~~~~~~~~~~~~~~~~~~~~~~Parameters
save = True
show = True

usepvalhack = False

fitfunc = twobody1

mapder = ''
zmin = 0
zmax = 0

mapcolor = 'Greys'#'plasma'
linecolor = 'brg'#'plasma'

gate = -6.4
source = [0,2.5]
modnumber = 2   #number of curves to plot
# source = [[.5, 1], [.5, 1], [.5, 1]]
powerdisplayvalue = .15 #power display value
smoothmap = True    #does the display map get smoothed?
smoothmapvalue = 0

smooth = False
smoothval = 0
smoothpower = False
smoothpval = 0
filtery = False
filterparams = (3, .18)
zeroset = True  #shifts x so x[0] is zero + poffset
sdzero = [-1, 3]
gatezero = [-4,-2]
fullzeroset = False #sets universal zero based on sd and gate location
pzeroset = True #shifts y so y[0] is zero + poffset
addzero = False #adds 0,0 as a starting value
falsezeroset = False    #sets the lowest power map to entirely zero
doffset = 0.001#1.994e-8
poffset = 0.001
lowpindex = 0
lowcindex = 2
trimval = .25
errorval = 1
p0 = [100,.01]
trackp0 = False

dabs = False    #absolute value
dnorm = False   #norms curves to 1

trimaxes = False    #plotting stuff
alim = (-1, 5000)
alim = (1e-7, 1e11)
alim = (-1, 1e11)
alim = (0)
alog = False
showfits = True
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

if usepvalhack:
    pval = getpvalhack(pval)
else:
    pval = []
    for i in range(plen):
        pval.append(np.mean(powdata[:,:,i]))
    pval = np.asanyarray(pval)

plen = pval.size

if fullzeroset:
    sdzeroindex = np.searchsorted(xval, sdzero)
    gatezeroindex = np.searchsorted(yval, gatezero)
    dataoffset = np.mean(data[gatezeroindex[0]:gatezeroindex[1],sdzeroindex[0]:sdzeroindex[1],0])
    data = data - dataoffset

powerdisplayindex = np.searchsorted(pval, powerdisplayvalue)
mapdata = data[:,:,powerdisplayindex]

if smoothmap:
    mapdata = applyGauss(mapdata, smoothmapvalue)

if mapder == 'x':
    derstep = np.abs(xval[0] - xval[1])
    ydelta, mapdata = np.gradient(mapdata[:,:], derstep)

if addzero:
    val0 = np.zeros(data[:,:,0].shape)
    pval = np.insert(pval, 0,0)
    data = np.insert(data, 0, val0,2)
    powdata = np.insert(powdata, 0, val0,2)
    plen += 1

if falsezeroset:
    data[:,:,0] = np.zeros(data[:,:,0].shape)

if smooth:
    for i in range(plen):
        data[:,:,i] = applyGauss(data[:,:,i], smoothval)

#~~~~~~~~~~~~~~~~~~~~~~Plot Generation

axs, axcbs, axscbs, figs = makeaxes(5, cbar = True, scalebar=True, scale = .7)
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


#plots mapdata
cmap, cnorm = color_memap(mapcolor, mapdata[:,:] , dmin = zmin, dmax = zmax)
xt = get_extent(xval, yval, inverty=True)

ax[0].imshow(mapdata[:,:] , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

# for i in range(len(gate)):
#     ax[0].axhline(gate[i], linestyle = ":", c = 'k', linewidth = 1)
ptrim = np.searchsorted(pval, trimval)
etrim = np.searchsorted(pval, errorval)

#tries all the fits

sdx = []
alpha = []
beta = []
c = []
sdc = []
error = []

if filtery:
    b, a = butter(filterparams[0], filterparams[1])

if yval[1] > yval[0]:
    gateindex = np.searchsorted(yval, gate)
else:
    gateindex = np.searchsorted(-yval, -gate)
if xval[1] > xval[0]:
    lowsdindex = np.searchsorted(xval, source[0])
    highsdindex = np.searchsorted(xval, source[1])
    sdcount = highsdindex - lowsdindex
else:
    lowsdindex = np.searchsorted(-xval, -source[0])
    highsdindex = np.searchsorted(-xval, -source[1])
    sdcount = lowsdindex - highsdindex 
clr = color_meline(linecolor, sdcount)
for p in range(sdcount):
    if p % modnumber == 0:
        y = data[gateindex, lowsdindex + p, :]
        x = powdata[gateindex, lowsdindex + p, :]

        if pzeroset:
            x = x - np.min(x) + poffset
            # print(offset)
        if zeroset:
            y = y - y[0] - doffset
        else:
            y = y + doffset

        if dabs:
            y = np.abs(y)

        if dnorm:
            scaling = 1 / np.abs(np.mean(y[ptrim-5:ptrim]))
            y = y * scaling

        # if smoothpower:
        #     y = applyGauss(y, smoothpval)

        if filtery:
            y = filtfilt(b, a, y)
        ax[1].plot(x[:ptrim], y[:ptrim], c = clr[p], linestyle = "-", linewidth = 1.4)
        ax[1].plot(x[ptrim:], y[ptrim:], c = clr[p], linestyle = ":", linewidth = 1.4)
        ax[0].scatter(xval[lowsdindex + p], gate, s = 2, c = clr[p] )
        try:
            if trackp0:
                par, pcov = curve_fit(fitfunc, x[lowpindex:ptrim], y[lowpindex:ptrim], p0 = p0, maxfev = 3200)
                p0 = par
            else:
                par, pcov = curve_fit(fitfunc, x[lowpindex:ptrim], y[lowpindex:ptrim], maxfev = 3200)
                if par[0] < 0:
                    par[0] = 0
            # error2b = 0
            # for qq in range(1, etrim, 1):
            #     error2b += np.abs( ( y[qq] - fitfunc(x[qq], *par ) / y[qq] ) )
            error2b = (pcov[0,0])
            parl, pcovl = curve_fit(line, x[:lowcindex], y[:lowcindex], p0 = [1,0], maxfev = 3200)
            # vary = 3*np.std(y[:ptrim])
            # ry = np.abs(np.mean(y[ptrim-4:ptrim]) - np.mean(y[:4]))
            vary = 1.2*np.abs(np.max(y[y.size//2:]) - np.min(y[y.size//2:]))
            vary = 2*np.abs(np.max(y[ptrim-10:ptrim]) - np.min(y[ptrim-10:ptrim]))
            ry = np.abs(np.mean(y[ptrim-3:ptrim]) - np.mean(y[:3]))
            if ry < vary:
                ax[1].plot(x[:ptrim], y[:ptrim], c = 'k', linestyle = "-", linewidth = 1.4)
                par[0] = 0
                par[1] = 0
            if showfits:
                ax[1].plot(x[:ptrim], fitfunc(x[:ptrim], *par), c = clr[p], linewidth = 1)
                # ax[1].plot(x[:ptrim], line(x[:ptrim], *parl), c = clr[p], linewidth = 1)
            ax[2].scatter(xval[lowsdindex + p], par[0], c = [clr[p]], s = 5)
            ax[3].scatter(xval[lowsdindex + p], par[1], c = [clr[p]], s = 5)
            # ax[4].scatter(xval[lowsdindex + p], parl[0], c = [clr[p]], s = 5)
            ax[4].scatter(xval[lowsdindex + p], vary, c = [clr[p]], s = 5)
            ax[4].scatter(xval[lowsdindex + p], ry, c = 'k', s = 5)
            # ax[4].scatter(xval[lowsdindex + p], error2b, c = [clr[p]], s = 5)
            sdx.append(xval[lowsdindex +p])
            alpha.append(par[0])
            error.append(error2b)
            beta.append(par[1])
            c.append(parl[0])
            # print(par[0])
        except RuntimeError:
            sdx.append(xval[lowsdindex +p])
            alpha.append(0)
            error.append(np.inf)
            beta.append(0)
            c.append(0)

        # try:
            
            
        #     sdc.append(xval[p])
        # except RuntimeError:
        #     None
# ax[1].plot(x, twobody(x, *p0), c = 'r', linewidth = 4)
# ax[1].set_xlim(np.min(pval[:ptrim]), np.max(pval[:ptrim]))
ax[1].set_title("gate = " + str(gate) + " V")
ax[1].set_xlabel("power (mW)")
ax[1].set_ylabel("Photocurrent (nA)")
barnorm = mpl.colors.Normalize(vmin = np.max(xval), vmax = np.min(xval))
mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(linecolor), norm = barnorm, orientation='horizontal')
    
sdx = np.asanyarray(sdx)
alpha = np.asanyarray(alpha)

# ax[2].scatter(sdx, alpha, s = 5, c = 'k')
# ax[3].scatter(sdx, beta, s = 5, c = 'k')
# ax[4].scatter(sdx, c, s = 5, c = 'k')

for i in range(alpha.size):
    print(str(sdx[i]) + "  :  " + str(alpha[i])+ "  :  " + str(error[i]))

if trimaxes:
    ax[1].set_xlim((pval[0], trimval))

if alog:
    ax[2].set_yscale('log')

if alim != [0]:
    ax[2].set_ylim(alim)

for i in range(3):
    ax[2+i].axhline(0, linestyle = ":", linewidth = 1, c = "k")
    ax[2+i].set_xlabel("Vsd (V)")

ax[3].set_title("Beta / A")
ax[2].set_title("Alpha / Gamma")
ax[4].set_title("Error **2")
ax[0].set_xlabel("Vsd (V)")
ax[0].set_xlabel("Power (mW)")

if save:
    savename = name + "_gate_trimpowercuts"
    pl.savefig("C:/Jed/Figures/" + savename + ".png")

if show:
    pl.show()