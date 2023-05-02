'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
from scipy.signal import butter, filtfilt
from scipy import fftpack
from datetime import date

from plot_tools import *
from analysis_tools import *
import loader

# for i in range(9, 31):

# name = "2021_03_25_51_powercube_twobody_" + str(i)
# name = "2021_03_25_51_powercube_twobody_" + str(ptrim) + "_" + str(smooth)
name2 = "2021_03_25_51_powercube_twobody_26_0.5"
name = "2021_03_25_51_powercube_twobody_28_0.5"
name1 = "2021_03_25_51_powercube_twobody_29_0.5"
name3 = "2021_03_25_51_powercube_twobody_27_0.5"
# name3 = "2021_03_25_51_powercube_twobody_26"
path = "C:/QMO/Built"

alpha, beta, cmap, alphaerror, betaerror, cerror, noise, yrange, xval, yval, rfi, info = loader.load_twobodyfit(path, name)
alpha1, beta1, cmap1, alphaerror1, betaerror1, cerror1, noise1, yrange1, xval1, yval1, rfi1, info1 = loader.load_twobodyfit(path, name1)
alpha2, beta1, cmap1, alphaerror1, betaerror1, cerror1, noise1, yrange1, xval1, yval1, rfi1, info1 = loader.load_twobodyfit(path, name2)
alpha3, beta1, cmap1, alphaerror1, betaerror1, cerror1, noise1, yrange1, xval1, yval1, rfi1, info1 = loader.load_twobodyfit(path, name3)
# alpha = (alpha + alpha1 + alpha2) * 0.3333
alpha = (alpha + alpha1) * 0.5
noise = (noise + noise1) / 2
yrange = (yrange + yrange1) / 2
alphaerror = (alphaerror + alphaerror1) / 2
betaerror = (betaerror + betaerror1) / 2
# alpha = (alpha + alpha1 + alpha2 + alpha3) * 0.25
# alpha, beta, cmap = alphaerror, betaerror, cerror
# noise = noise * 2
error = np.sqrt(alphaerror * alphaerror + betaerror * betaerror)
cmap = error
info = make_dict(info)

xlabel = "V_{sd} (V)"
ylabel = "V_{gate} (V)"
# xlabel = ""
# ylabel = ""

savemod = ""
save = True
show = True

mapcolor_name = ['magma', 'plasma', 'viridis', 'seismic', 'plasma'] #Color of map
colormap_name = 'viridis'

smooth = [.7, 0, 0]
upperlimit = [10000000000, 0, 0]
errortrim = 0
hotpixels = False
hotpixelstep = 5
lowpass = False
order = 3
cutoff = 2/28 # cutoff freq / .5 * sample rate

noisetrim = 2
noisetrimsdcutoff = .9

xslicelimits = [[],[],[],[],[]]
yslicelimits = [[0, 2000],[],[],[],[]]
ymaplimits = []
xmaplimits = []
zlimit = [[0, 2300],[],[],[],[]] 

hlines = []
vlines = []

slices = 50
gateCuts = []
sourceCuts = []

gatemin = -60
gatemax = 50
sdmin = -50
sdmax = 50

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Data manipulation

xlen = xval.size
ylen = yval.size

if errortrim > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(error[p,i]) > errortrim:
                alpha[p,i] = 0
            if alpha[p,i] < 0:
                alpha[p,i] = 0
            if alpha[p,i] == np.inf or alpha[p,i] == -np.inf:
                alpha[p,i] = 0
            
if noisetrim > 0:
    for i in range(xlen):
        for p in range(ylen):
            if xval[i] < noisetrimsdcutoff:
                if noisetrim * noise[p,i] > yrange[p,i]:
                    alpha[p,i] = 0

if hotpixels:
    for n in range(ylen):
        for m in range(xlen-1):
            # print( str(np.abs(alpha[n,xlen-m-2])) + "  " + str(np.abs(alpha[n,xlen-m-2])))
            if np.abs(alpha[n,xlen-m-2]) > hotpixelstep * np.abs(alpha[n,xlen-m-1]):
                alpha[n,xlen-m-2] = alpha[n,xlen-m-1]
    # for n in range(ylen):
    #     for m in range(xlen-1):
    #         # print( str(np.abs(alpha[n,xlen-m-2])) + "  " + str(np.abs(alpha[n,xlen-m-2])))
    #         if np.abs(alpha[n,m+1]) > hotpixelstep * np.abs(alpha[n,m]):
    #             alpha[n,m+1] = alpha[n,m]


# pl.figure()
# clr = color_meline('plasma', ylen)
# y = alpha[10,:]
# signalnoise = fftpack.fft(y)
# signalamp = xlen * np.abs(signalnoise)
# signalfreq = np.abs(fftpack.fftfreq(xlen, 2 / xlen))
# pl.plot(signalfreq, signalamp, linewidth = .5, c = clr[50])
# pl.show()

if lowpass:
    b, a = butter(order, cutoff, analog=False)
    for i in range(ylen):
        dy = alpha[i,:]            
        y = filtfilt(b, a, dy)
        alpha[i,:] = y[:]

if upperlimit[0] > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(alpha[p,i]) > upperlimit[0]:
                alpha[p,i] = 0
if upperlimit[1] > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(beta[p,i]) > upperlimit[1]:
                beta[p,i] = 0
if upperlimit[2] > 0:
    for i in range(xlen):
        for p in range(ylen):
            if np.abs(cmap[p,i]) > upperlimit[2]:
                cmap[p,i] = 0

if gateCuts == []:
    gateCuts = np.linspace(np.min(yval), np.max(yval), slices)
if sourceCuts == []:
    sourceCuts = np.linspace(np.min(xval), np.max(xval), slices)

if smooth[0] > 0:
    alpha = applyGauss(alpha, smooth[0])
if smooth[0] > 0:
    beta = applyGauss(beta, smooth[1])
if smooth[0] > 0:
    cmap = applyGauss(cmap, smooth[2])

for i in range(len(zlimit)):
    if zlimit[i] == []:
        zlimit[i] = [0,0]

data = [alpha, beta, cmap, 2*noise - yrange, yrange]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Figure Generation
print("Building Figures")
sets = 5
axs, figs = make_grid(3,sets)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Plot Population
print("Populating Plots")

xt = get_extent(xval, yval, inverty= True)
for i in range(sets):
    cmap, cnorm = color_memap(mapcolor_name[i], data[i] , dmin = zlimit[i][0], dmax = zlimit[i][1])
    ax[3*i].imshow(data[i], cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')


if gateCuts != []:
    clrgate = color_meline(colormap_name, len(gateCuts)  )
if sourceCuts != []:
    clrsd = color_meline(colormap_name, len(sourceCuts)  )

for q in range(sets):
    for i in range(len(gateCuts)):
        if gateCuts[i] > gatemin and gateCuts[i] < gatemax:
            gateindex = np.searchsorted(yval, gateCuts[i])
            x = xval
            y = data[q][gateindex, :]
            ax[3*q+1].plot(x, y, linewidth = 1.5, c = clrgate[i])
    for i in range(len(sourceCuts)):
        if sourceCuts[i] > sdmin and sourceCuts[i] < sdmax:
            sdindex = np.searchsorted(xval, sourceCuts[i])
            x = yval
            y = data[q][:, sdindex]
            ax[3*q+2].plot(x, y, linewidth = 1.5, c = clrsd[i])
    
    ax[3*q+1].set_xlim(xval[0], xval[-1])
    ax[3*q+2].set_xlim(yval[0], yval[-1])

    if yslicelimits[q] != []:
        ax[3*q+1].set_ylim(yslicelimits[q][0], yslicelimits[q][1])
        ax[3*q+2].set_ylim(yslicelimits[q][0], yslicelimits[q][1])

    ax[3*q+1].set_title("Horizontal Cuts")
    ax[3*q+2].set_title("Vertical Cuts")

ax[0].set_title("Alpha")
ax[3].set_title("Beta")
ax[6].set_title("CMap")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Save/Show
print("Saving")
if save:
    today = date.today()
    savename = name + "__" + str(today)
    pl.savefig("C:/QMO/generated plots/" + savemod + "/" + savename +".png")
if show:
    pl.show()
