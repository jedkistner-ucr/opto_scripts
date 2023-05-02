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
power = 400
alpha = [0.000001, .15]
alpha = [0.000001, .06]
# alpha = [0.000001, .15]
tau = .001

alpha = np.linspace(alpha[0], alpha[1], 100)
arate = alpha / tau

axs, figs = customPlot(5, 5)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))


y = twobody(alpha, bvalue, power)*10e9
ymax = np.max(y)
ax[0].plot(arate, y, color = 'k')

ax[0].set_xlabel("A (a/t)")
ax[0].set_ylabel("Photocurrent at 400mW (nA)")
ax[0].set_xlim(-2, np.max(arate))
ax[0].axvline(0, linestyle = ":", c = 'k')
ax[0].set_ylim(0, ymax * 1.1)

alpha_ = [0.000001, .011]
tau_ = .001
curvenumber = 10
alpha_ = np.linspace(alpha_[0], alpha_[1], curvenumber)
arate_ = alpha_/tau_
clr = color_meline('plasma', curvenumber)
for i in range(alpha_.size):
    ax[0].axvline(arate_[i], linestyle = ":", c = clr[i])

pl.savefig("C:/QMO/generated plots/electron_hole/" + "current_sensitivity" +".png")

pl.show()

