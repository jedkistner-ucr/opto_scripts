import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import random
import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *

#energy pigment a
u_a = 1.1
p_a = .5

#energy pigment b
u_b = 0.9
p_b = .5

#machine energy and timestep
u_m = 1
dt = 1

rounds = 20000
ext = 0

# energy = 0
eval = []
tval = []

# x = np.asarray(range(rounds))
# ext = .01 * np.sin(x/50)
pigs = 10

for i in range(rounds):
    energy = 0
    if i % 20 == 0:
        ext = ((random.random()) * 0.05)
        signflip = random.random()
        if signflip > 0.5:
            ext = -ext
    for p in range(pigs):
        r = random.random()
        #add energy
        if r < p_a:
            energy += (u_a + ext) * dt
        elif r > p_a and r < (p_a + p_b):
            energy += (u_b + ext) * dt
    
    #subtract energy
    # if energy > (10* u_m * dt):
    #     energy = energy - (10 * u_m * dt)

    tval.append(i)
    eval.append(energy)

tval = np.asarray(tval)
eval = np.asarray(eval)

eval = eval / pigs

# axs, figs = customPlot(4, 7, sidebar = False)
# ax = []
# fig = pl.figure(num = 0, figsize=figs)
# for a in axs:
#     ax.append(fig.add_axes(a))

# ax[0].axvline(1, c = 'k', linewidth = 1, linestyle = ":")
# ax[0].plot(eval, tval, c = 'b', linewidth = 1)
# ax[0].set_ylim(tval[-1], tval[0])
# xmargin = np.max(np.abs(eval - 1)) * 1.1
# ax[0].set_xlim(-xmargin+1, xmargin+1)
# ax[0].yaxis.tick_right()
# ax[0].xaxis.tick_top()
# ax[0].set_xticks([.9, 1, 1.1])

# bins = np.histogram(eval, bins = 21)
# pl.scatter(bins[1][1:], bins[0])

axs, figs = customPlot(4, 3, sidebar = False)
ax1 = []
fig1 = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax1.append(fig1.add_axes(a))

ax1[0].yaxis.tick_right()
# ax1[0].hist(bins[0], 21)
bins_ = np.linspace(.85, 1.15, 22)
n, bins, patches = ax1[0].hist(eval, bins_)
# ax1[0].scatter(bins[1][1:], bins[0])
clr = color_meline('plasma', 21)
pq = 0
for p in patches:
    if pq == 10:
        pl.setp(p, 'facecolor', 'orange')
    elif pq < 10:
        pl.setp(p, 'facecolor', 'blue')
    else:
        pl.setp(p, 'facecolor', 'red')
    pl.setp(p, 'edgecolor', 'k')
    pq += 1

xmargin = np.max(np.abs(eval - 1)) * 1.1
# ax1[0].set_xlim(-xmargin+1, xmargin+1)
ax1[0].set_xlim((.85, 1.15))
ax1[0].set_xticks([1])

# ax1[0].yaxis.tick_right()
savename = "histv3"
pl.savefig("C:/QMO/photosynthesis/" + savename + ".png")
pl.show()
    
