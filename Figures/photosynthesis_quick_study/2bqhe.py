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

rounds = 5000
ext = 0

executions = 50
pigs = 10

trim = 1000

for x in range(executions):

    eval = []
    tval = []

    for i in range(rounds):
        energy = 0
        if i % 50 == 0:
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

    axs, figs = photosyn_plot()
    ax = []
    fig = pl.figure(figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))


    ax[0].plot(eval[:trim], tval[:trim], c = 'b', linewidth = 1)
    ax[0].set_ylim(tval[trim], tval[0])
    xmargin = np.max(np.abs(eval - 1)) * 1.1

    ax[0].set_xlim((.85, 1.15))
    ax[0].yaxis.tick_right()
    ax[0].xaxis.tick_top()
    ax[0].set_xticks([.9, 1, 1.1])

    ax[1].yaxis.tick_right()
    bins_ = np.linspace(.85, 1.15, 22)
    n, bins, patches = ax[1].hist(eval, bins_)

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

    ax[1].set_xlim((.85, 1.15))
    ax[1].set_xticks([1])

    savename = "fig1_2_" + str(x)
    pl.savefig("C:/QMO/photosynthesis/" + savename + ".png")
    if executions > 1:
        pl.cla()

if executions == 1:
    pl.show()
    
