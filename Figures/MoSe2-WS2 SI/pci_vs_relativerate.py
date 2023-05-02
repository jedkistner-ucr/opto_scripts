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

def dn_dt(alpha, tau, n):
    a = alpha * n * n
    b = tau * n
    return -a, -b

alpha = [0, .011]
tau = .001

alpha = np.linspace(alpha[0], alpha[1], 200)
arate = alpha / tau
pci_tot = np.zeros(arate.size)

initialpop = 100
timesteps = 300
deltat = 1
laserstep = 100

for p in range(arate.size):

    t = np.arange(timesteps)
    N = np.zeros(t.size)
    pci = np.zeros(t.size)
    eh = np.zeros(t.size)
    N[0] = initialpop

    for i in range(timesteps-1):
        a, b = dn_dt(alpha[p], tau, N[i])

        N[i+1] = ( (a+b) * deltat) + N[i]
        pci[i] = (-b * deltat)
        eh[i] =  (-a * deltat)
        # pci_tot[i+1] = pci_tot[i] + (-b * deltat)

        if N[i+1] < 0:
            N[i+1] = 0
            pci[i] = 0
            eh[i] =  0

        if i % laserstep == 0 and i > 0:
            N[i+1] = N[i+1] + initialpop

    pci_tot[p] = np.sum(pci)

axs, figs = customPlot(6, 4)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

# ax0 = ax[0].twinx()
# ax1 = ax[0].twinx()

ax[0].plot(arate, pci_tot, color = 'k')
# ax0.plot(t, pci, color = 'r')
# ax1.plot(t, pci_tot, color = 'r', linestyle = ":")
# ax0.plot(t, eh, color = 'b')

# ax1.set_ylim(0, 65)
# ax[0].set_ylim(0, 120)
# ax0.set_ylim(0, 120)

# ax[0].set_xlim(0,1)

# ax0.set_yticks([])

ax[0].set_xlabel("Relative Interaction Rate (unitless)")
ax[0].set_ylabel("Photocurrent (Arb. Amps)")
# ax1.set_ylabel("Total photocurrent produced (arb. Amps)")
# ax1.yaxis.label.set_color('red')
# ax1.tick_params(axis='y', colors='red')

# ax[0].set_title("alpha = %.1e  tau = %.1e" % (alpha, tau))

# pl.savefig("C:/QMO/generated plots/electron_hole/" + "www" +".png")

pl.show()

