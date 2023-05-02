import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import math
import csv

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *


sunpath = "C:/QMO/solar.csv"

#blackbody radiation equation. Takes an array of frequencies and returns the radiance across that spectrum
def blackbody(t, w):
    c =  3*10e8 #m/s
    h = 6.626070*10e-34 #j/hz
    kt = t * (1.380649*10e-23) #j

    # w = w*10e-8

    # radiance = ((2*h*c*c)/(w**5)) * ( 1/ (np.exp( (h*c)/(w*kt) ) - 1) )
    a = 2.0*h*c**2
    b = h*c/(w*kt)
    radiance = a/ ( (w**5) * (np.exp(b)-1) )
    # x = c/w 
    # radiance = ((2*h*x*x*x)/(c*c))*(1/(np.exp(h*x/kt)-1))
    return radiance

def lorentz(center, width, x):
    y = (width / 2*math.pi)*(1 / ( ((x - center)**2) + ((width/2)**2) ))
    return y

#frequency range and temperatures to calculate for
save = True
w = np.linspace(200, 800, 1200)
t = [5778]
colormap = 'gist_rainbow_r'
colormap1 = 'bwr_r'
colormap2 = 'bwr_r'
clr = ['k']

vl = [400e-9, 600e-9]
vl=[]

smallrange = [608*10e-9, 808*10e-9]
smallrange_ = [200*10e-9, 400*10e-9]

lc = [668*10e-9, 748*10e-9, 336*10e-9, 256*10e-9]
lc = [429, 459]
lw = [10, 10]
mod = 3
submod = -.2
clr_ = ['b', 'c', 'r', 'b']

xlimit = (2e-6, 8e-6)
ylimit = (0, 1.2)
w = w * 10e-9

xti = np.linspace(0, 1200*10e-9, 7)
xla_ = np.linspace(0, 1200, 7)
xla = []
for a in xla_:
    xla.append(int(a))

sundata = np.genfromtxt(sunpath, delimiter = ",")




#Plot generation
axs, figs = customPlot_topbar(5, 4)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
#color generation
if len(clr) == 0:
    clr = color_meline(colormap, len(t))

maxval = 0


rad = blackbody(t[0], w)
rad = rad / np.max(rad)
radgrad = np.gradient(rad)
minindex = np.argmin(radgrad)
maxindex = np.argmax(radgrad)

# w = sundata[:,0]
# rad = sundata[:,2]
ylimit = (0, np.max(rad) * 1.2)
zeroline = np.zeros(w.shape)

# ax[0].axvline(w[minindex], c = 'k', linewidth = 2, alpha = 0.5, linestyle = ":")
# ax[0].axvline(w[maxindex], c = 'k',  linewidth = 2, alpha = 0.5, linestyle = ":")
ax[0].plot(w, rad, linewidth = 2, c = clr[0], alpha = 1)
# for v in vl:
#     ax[0].axvline(v, c = 'k', linewidth = 2, alpha = 0.5, linestyle = ":")
# ax[0].axvline(w[maxindex], c = 'k',  linewidth = 2, alpha = 0.5, linestyle = ":")

lowrange = np.searchsorted(w, smallrange[0])
highrange = np.searchsorted(w, smallrange[1])
lowrange_ = np.searchsorted(w, smallrange_[0])
highrange_ = np.searchsorted(w, smallrange_[1])

# ax[0].plot(w[lowrange:highrange], rad[lowrange:highrange]+0e10, linewidth = 2.1, c = 'orange', alpha = 1)
# ax[0].plot(w[lowrange_:highrange_], rad[lowrange_:highrange_]+0e10, linewidth = 2.1, c = 'orange', alpha = 1)
# ax[0].plot(w, radgrad, linewidth = 2, c = clr[0], alpha = 0.7)


for i in range(len(lc)):
    y = (lorentz(lc[i], lw[i], w) * mod) + submod
    ax[0].plot(w, y, linewidth = 1.5, c = clr_[i], alpha = 1)

    # for p in range(y.size):
    #     if y[p] > rad[p]:
    #         y[p] = rad[p]
    # fillx = np.concatenate((w, w))
    # filly = np.concatenate((zeroline, y))
    # ax[0].fill(fillx, filly, c = clr_[i], alpha = 0.4)
ax[-3].spines.right.set_visible(False)
ax[-2].spines.left.set_visible(False)
ax[-1].spines.right.set_visible(False)
ax[-1].spines.left.set_visible(False)

cmap, cnorm = color_memap(colormap, np.asarray([-1,1]) , dmin = 0, dmax = 0)
mpl.colorbar.ColorbarBase(ax[-1], cmap = cmap, norm = cnorm, orientation='horizontal')

cmap, cnorm = color_memap(colormap1, np.asarray([-1,1]) , dmin = 0, dmax = 0)
mpl.colorbar.ColorbarBase(ax[-2], cmap = cmap, norm = cnorm, orientation='horizontal')
cmap, cnorm = color_memap(colormap2, np.asarray([-1,1]) , dmin = 0, dmax = 0)
mpl.colorbar.ColorbarBase(ax[-3], cmap = cmap, norm = cnorm, orientation='horizontal')
ax[-1].set_xticks([])
ax[-1].set_xlim(-.57, 1)
ax[-2].set_xlim(-1,0)
ax[-3].set_xlim(0,1)

ax[-2].set_xticks([])
ax[-2].set_yticks([])
ax[-3].set_xticks([])
ax[-3].set_yticks([])

ax[0].xaxis.tick_top()
ax[0].xaxis.set_label_position('top') 

ax[0].set_xticks([2e-6,4e-6, 6e-6, 8e-6])

ax[0].set_xlim(xlimit)
if ylimit != ():
    ax[0].set_ylim(ylimit)

if save:
    pl.savefig("C:/QMO/generated plots/photosyn/" + "blackbody_2" +".png")

pl.show()
