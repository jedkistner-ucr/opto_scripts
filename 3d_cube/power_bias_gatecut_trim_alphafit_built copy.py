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
name = '2021_03_25_51_powercube_butter_4_3'
name = '2021_03_25_51_powercube_butter_2_4'
path = "C:/Jed/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
# data = data * 1e9
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

gate = -6.5
source = [0,2.5]
modnumber = 1   #number of curves to plot
# source = [[.5, 1], [.5, 1], [.5, 1]]
powerdisplayvalue = .15 #power display value
smoothmap = True    #does the display map get smoothed?
smoothmapvalue = 0

smooth = True
smoothval = 1
smoothpower = False
smoothpval = 0
zeroset = True  #shifts x so x[0] is zero + poffset
pzeroset = True #shifts y so y[0] is zero + poffset
addzero = False #adds 0,0 as a starting value
falsezeroset = False    #sets the lowest power map to entirely zero
doffset = 0
poffset = 0
# trimval = .15
ptrimindex = 28

p0 = [100,.01]

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

if usepvalhack:
    pval = getpvalhack(pval)
else:
    pval = []
    for i in range(plen):
        pval.append(np.mean(powdata[:,:,i]))
    pval = np.asanyarray(pval)

plen = pval.size

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

ptrim = ptrimindex

#tries all the fits

sdx = []
alpha = []
beta = []
r = []

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
           
        if zeroset:
            y = y - y[0] - doffset

        else:
            y = y + doffset

        if dabs:
            y = np.abs(y)

        if dnorm:
            scaling = 1 / np.abs(np.mean(y[ptrim-5:ptrim]))
            y = y * scaling

        ax[1].plot(x[:ptrim], y[:ptrim], c = clr[p], linestyle = "-", linewidth = 1.4, alpha = .2)
        ax[1].plot(x[ptrim:], y[ptrim:], c = clr[p], linestyle = ":", linewidth = 1.4, alpha = .2)
        ax[0].scatter(xval[lowsdindex + p], gate, s = 2, c = clr[p] )
        try:
            par, pcov = curve_fit(fitfunc, x[:ptrim], y[:ptrim], maxfev = 3200)
            if par[0] < 0 or par[0] > 10000:
                par[0] = 0

            res_squared = np.sum((y-fitfunc(x,*par))**2)
            sum_squared = np.sum((y-np.mean(y))**2)
            r_val = 1 - (res_squared/sum_squared)
            
            if showfits:
                ax[1].plot(x[:ptrim], fitfunc(x[:ptrim], *par), c = clr[p], linewidth = 1)
                # ax[1].plot(x[:ptrim], line(x[:ptrim], *parl), c = clr[p], linewidth = 1)
            ax[2].scatter(xval[lowsdindex + p], par[0], c = [clr[p]], s = 5)
            ax[3].scatter(xval[lowsdindex + p], par[1], c = [clr[p]], s = 5)
            ax[4].scatter(xval[lowsdindex + p], r_val, c = 'k', s = 5)

            sdx.append(xval[lowsdindex +p])
            alpha.append(par[0])
            r.append(r_val)
            beta.append(par[1])


        except RuntimeError:
            sdx.append(xval[lowsdindex +p])
            alpha.append(0)
            r.append(0)
            beta.append(0)



# ax[1].plot(x, twobody(x, *p0), c = 'r', linewidth = 4)
# ax[1].set_xlim(np.min(pval[:ptrim]), np.max(pval[:ptrim]))
ax[1].set_title("gate = " + str(gate) + " V")
ax[1].set_xlabel("power (mW)")
ax[1].set_ylabel("Photocurrent (nA)")
barnorm = mpl.colors.Normalize(vmin = np.max(xval), vmax = np.min(xval))
mpl.colorbar.ColorbarBase(axsc[0], cmap = pl.get_cmap(linecolor), norm = barnorm, orientation='horizontal')
    
sdx = np.asanyarray(sdx)
alpha = np.asanyarray(alpha)

# if trimaxes:
#     ax[1].set_xlim((pval[0], pval[ptrim]))
# else:
#     ax[1].set_xlim((pval[0], pval[-1]))

# ax[1].set_ylim(0)

ax[4].set_ylim(0, 1.05)

if alog:
    ax[2].set_yscale('log')

if alim != [0]:
    ax[2].set_ylim(alim)

for i in range(3):
    ax[2+i].axhline(0, linestyle = ":", linewidth = 1, c = "k")
    ax[2+i].set_xlabel("Vsd (V)")

ax[3].set_title("Beta / A")
ax[2].set_title("Alpha / Gamma")
ax[4].set_title("R value")
ax[0].set_xlabel("Vsd (V)")
ax[0].set_xlabel("Power (mW)")


if save:
    savename = name + "_gate_trimpowercuts"
    pl.savefig("C:/Jed/Figures/" + savename + ".png")

if show:
    pl.show()