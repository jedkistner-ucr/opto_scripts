'''
Models the excitation equation with time
'''
from time import time
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *

def twobody(A, B, P):
    f = 76e6
    deltat = .003
    q = 1.60217663e-19

    i = q * (B/A) * np.log(1+A*P) * f * deltat
    return i


bvalue = 1
prange = [0, 400]
alpha = [0.000001, .011]
tau = .001
curvenumber = 10

power = np.linspace(prange[0], prange[1], 100)

alpha = np.linspace(alpha[0], alpha[1], curvenumber)
clr = color_meline('plasma', curvenumber)
arate = alpha / tau

axs, figs = customPlot(5, 5)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

ymax = 0

for i in range(curvenumber):
    y = twobody(alpha[i], bvalue, power)*10e9
    if np.max(y) > ymax:
        ymax = np.max(y)
    ax[0].plot(power, y, color = clr[i])

ax[0].set_xlabel("Power (mW)")
ax[0].set_ylabel("Photocurrent (nA)")
ax[0].set_xlim(prange[0], prange[1])
ax[0].set_ylim(0, ymax * 1.1)

pl.savefig("C:/QMO/generated plots/electron_hole/" + "current_power_curves" +".png")

pl.show()

