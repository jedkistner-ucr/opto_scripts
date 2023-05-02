import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl

from datetime import date

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *

colormap = 'gist_rainbow'

axs, figs = customPlot_noedges(10, 1)
ax = []
fig = pl.figure(num = 0, figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

data = np.linspace(0, 1, 10)
ax[0].set_xticks([])

cmap, cnorm = color_memap(colormap, data , dmin = 0, dmax = 0)

mpl.colorbar.ColorbarBase(ax[0], cmap = cmap, norm = cnorm, orientation='horizontal')
ax[0].set_xticks([])

savename = "rainbow_bar"
pl.savefig("C:/QMO/photosynthesis/" + savename + ".png")

pl.show()