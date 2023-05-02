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
from curves import line, twobody, twobody1, simplepower

# def run(trim):

#~~~~~~~~~~~~~~~~~~~~~~Loading

name = "2022_02_03_12_powercube"
name = "2021_03_25_51_powercube"
path = "C:/QMO/Built"

xval, yval, pval, data, rf, powdata, info_ = loader.load_cube(path, name)
data = data * 1e9

saveadd = ""

# alphaerrorsave = "C:/QMO/Built/" + name + "_errormap0"+ saveadd
# betaerrorsave = "C:/QMO/Built/" + name + "_errormap1"+ saveadd
# betasave = "C:/QMO/Built/" + name + "_betamap"+ saveadd
# c_onesave = "C:/QMO/Built/" + name + "_cmap"+ saveadd

#~~~~~~~~~~~~~~~~~~~~~~Parameters
save = True

fitfunc = twobody

smooth = 0
smoothpower = 0


fullzeroset = False
sdzero = [-1, 3]
gatezero = [-4,-2]

zeroset_x = True
zeroset_y = True
addzero = False
falsezeroset = False

xoffset = 0
yoffset = 0
lowpindex = 0
lowcindex = 4
trimval = .23
errorval = 1
p0 = [10,.01]
trackp0 = False

dabs = False
dnorm = False

#~~~~~~~~~~~~~~~~~~~~~~Manipulation

xlen = xval.size
ylen = yval.size
plen = pval.size

if fullzeroset:
    sdzeroindex = np.searchsorted(xval, sdzero)
    gatezeroindex = np.searchsorted(yval, gatezero)
    dataoffset = np.mean(data[gatezeroindex[0]:gatezeroindex[1],sdzeroindex[0]:sdzeroindex[1],0])
    data = data - dataoffset

#Get noise data before smoothig anything
ptrim = np.searchsorted(pval, trimval)
noise = np.zeros(data[:,:,0].shape)
yrange = np.zeros(data[:,:,0].shape)
for i in range(ylen):
    for p in range(xlen):
        y = data[i, p, :ptrim]
        noise[i,p] = np.abs(np.max(y[-10:]) - np.min(y[-10:]))
        yrange[i,p] = np.abs(np.mean(y[-3:]) - np.mean(y[:3]))

if smooth > 0:
    for i in range(plen):
        data[:,:,i] = applyGauss(data[:,:,i], smooth)
#~~~~~~~~~~~~~~~~~~~~~~Plot Generation

alpha = np.zeros(data[:,:,0].shape)
beta = np.zeros(data[:,:,0].shape)
alphaerror = np.zeros(data[:,:,0].shape)
betaerror = np.zeros(data[:,:,0].shape)
cmap = np.zeros(data[:,:,0].shape)
cerror = np.zeros(data[:,:,0].shape)

blank = np.zeros(data[:,:,0].shape)

print("Starting fitting %i rows of points" % (ylen))

failedfits = 0

# print(powdata[0, 0, :])

for i in range(ylen):
    for p in range(xlen):
        y = data[i, p, :ptrim]
        x = powdata[i, p, :ptrim]

        if dabs:
            y = np.abs(y)

        if zeroset_x:
            x = x - x[0]
        if zeroset_y:
            y = y - y[0]
        
        x = x + xoffset
        y = y + yoffset

        if smoothpower > 0:
            y = applyGauss(y, smoothpower)

        if dnorm:
            scaling = 1 / np.abs(np.mean(y[ptrim-5:ptrim]))
            y = y * scaling

        try:
            if trackp0:
                if p == 0:
                    p0 = p0_
                par, pcov = curve_fit(fitfunc, x, y, p0 = p0, maxfev = 3200)
                p0 = par
                print(par[0])
            else:
                par, pcov = curve_fit(fitfunc, x, y, p0 = p0, maxfev = 3200)
            

            alpha[i,p] = par[0]
            beta[i,p] = par[1]
            alphaerror[i,p] = np.sqrt(np.diag(pcov))[0]
            betaerror[i,p] = np.sqrt(np.diag(pcov))[1]
        except RuntimeError:
            alpha[i,p] = 0
            beta[i,p] = 0
            alphaerror[i,p] = -1
            betaerror[i,p] = -1
            failedfits += 1

        try:
            parl, pcovl = curve_fit(line, x[:lowcindex], y[:lowcindex], p0 = [1,0], maxfev = 3200)
            cmap[i,p] = parl[0]
            cerror[i,p] = np.sqrt(np.diag(pcovl))[0]
        except:
            cmap[i,p] = 0
            cerror[i,p] = -1


    print(" " + str(i))

print("%i fits failed" % (failedfits))
print("  ~~ saving ~~  ")

if save:
    savename = "C:/QMO/Built/" + name + "_" + str(fitfunc).split(" ")[1] + saveadd
    np.savez(savename, alpha = alpha, beta = beta, cmap = cmap, alphaerror = alphaerror, betaerror = betaerror, cerror = cerror, noise = noise, yrange = yrange, xval = xval, yval = yval, rfi = rf[0,0,-1], info = info_)
    # np.savez(alphaerrorsave, d = alphaerror, xval = xval, yval = yval, rfi = blank, info = inf)
    # np.savez(betaerrorsave, d = betaerror, xval = xval, yval = yval, rfi = blank, info = inf)
    # np.savez(betasave, d = beta, xval = xval, yval = yval, rfi = blank, info = inf)
    # np.savez(c_onesave, d = cmap, xval = xval, yval = yval, rfi = blank, info = inf)

