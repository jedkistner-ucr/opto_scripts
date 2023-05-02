import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import math

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *


#blackbody radiation equation. Takes an array of frequencies and returns the radiance across that spectrum
def blackbody(t, w):
    c = 3*10e8 #m/s
    h = 6.626070*10e-34 #j/hz
    kt = t * (1.380649*10e-23) #j

    radiance = ((2*h*c*c)/(w**5)) * ( 1/ (np.exp( (h*c)/(w*kt) ) - 1) )
    return radiance

def lorentz(center, width, x):
    y = (width / 2*math.pi)*(1 / ( ((x - center)**2) + ((width/2)**2) ))
    return y

#frequency range and temperatures to calculate for
save = True
w = np.linspace(160, 1000, 1200)
t = [5778]
colormap = 'plasma'
clr = ['k']

smallrange = [608*10e-9, 808*10e-9]
smallrange_ = [200*10e-9, 400*10e-9]

lc = [668*10e-9, 748*10e-9, 336*10e-9, 256*10e-9]
lw = [40*10e-9, 40*10e-9, 40*10e-9, 40*10e-9]
mod = 2.2e4
submod = -.8e11
clr_ = ['r', 'b', 'r', 'b']

ylimit = (0, 3e11)
w = w * 10e-9

xti = np.linspace(0, 1200*10e-9, 7)
xla_ = np.linspace(0, 1200, 7)
xla = []
for a in xla_:
    xla.append(int(a))

zeroline = np.zeros(w.shape)

#Plot generation
axs, figs = customPlot(5, 5, sidebar = False)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))
#color generation
if len(clr) == 0:
    clr = color_meline(colormap, len(t))

maxval = 0


rad = blackbody(t[0], w)
radgrad = np.gradient(rad)
minindex = np.argmin(radgrad)
maxindex = np.argmax(radgrad)
# ax[0].axvline(w[minindex], c = 'k', linewidth = 2, alpha = 0.5, linestyle = ":")
# ax[0].axvline(w[maxindex], c = 'k',  linewidth = 2, alpha = 0.5, linestyle = ":")
ax[0].plot(w, rad, linewidth = 2, c = clr[0], alpha = 1)

lowrange = np.searchsorted(w, smallrange[0])
highrange = np.searchsorted(w, smallrange[1])
lowrange_ = np.searchsorted(w, smallrange_[0])
highrange_ = np.searchsorted(w, smallrange_[1])

# ax[0].plot(w[lowrange:highrange], rad[lowrange:highrange]+0e10, linewidth = 2.1, c = 'orange', alpha = 1)
# ax[0].plot(w[lowrange_:highrange_], rad[lowrange_:highrange_]+0e10, linewidth = 2.1, c = 'orange', alpha = 1)
# ax[0].plot(w, radgrad, linewidth = 2, c = clr[0], alpha = 0.7)


# for i in range(len(lc)):
#     y = (lorentz(lc[i], lw[i], w) * mod) + submod
#     ax[0].plot(w, y, linewidth = 1.5, c = clr_[i], alpha = 1)

#     for p in range(y.size):
#         if y[p] > rad[p]:
#             y[p] = rad[p]
#     fillx = np.concatenate((w, w))
#     filly = np.concatenate((zeroline, y))
#     ax[0].fill(fillx, filly, c = clr_[i], alpha = 0.4)

ax[0].set_xticks(xti, xla)

ax[0].set_xlim(w[0], w[-1])
if ylimit != ():
    ax[0].set_ylim(ylimit)

if save:
    pl.savefig("C:/QMO/generated plots/photosyn/" + "blackbody" +".png")

pl.show()
