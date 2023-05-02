'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

# import matplotlib as mpl
# import matplotlib.pyplot as pl
# import matplotlib.colors as colors
# from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map
from curves import simplepower, twobody1, exponential

def main(savemod = 0, plotgline = 0):

    name = 'A02_-4gate_sd'
    # name = "d02_deltaV"

    path = "C:/QMO/Built"


    xval, yval, data = loader.load_simplemap(path, name)
    # data = np.transpose(data)
    # yval = np.flip(yval)
    # data = np.flip(data, 1)
    # yval = 1240/yval


    xlabel = ""#"$V_{TG} = V_{BG}$  (V)"
    ylabel = ""#"$\it{I}$  (nA)"
    ylabel1 = ""#"$\it{I}$  (nA)"

    save = True
    show = True

    savemod = 1

    #Parameters
    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices
    contours = False

    colormap_name = "viridis_r"

    #How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
    slices = 300
    # gatecuts = np.linspace(-5, -6, 10)
    # gatecuts = np.linspace(0, 5, 8)
    # gatecuts = np.asarray([-1, -2, -3, -3.5, -4, -4.5, -5])
    # gatecuts = [1.3]
    # gatecuts = [2.5, -.5]
    # plotgline = [1, 1, 1, 1, 1, 1, 1]
    # plotgline = np.ones(gatecuts.shape)
    # plotgline = plotgline
    # speccolor = color_meline(colormap_name, gatecuts.size  )
    # speccolor = ["k", 'r', 'g']
    # style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
    alpha = 1
    width = 5
    clr = "k"

    gatecuts1 = [-5.13]
    plotgline1 = [1]
    speccolor1 = ["k"]
    style1 = ["-"]
    alpha1 = [1]
    width1 = [1]

    scaling = 1
    scaling1 = 1
    offset = -.15
    offset1 = 0

    direrction = "x"

    hlines = [1]
    vlines = [0]

    #simple parameters -- in the otf version these pretty much stay the way they are
    interpolate = False
    newbox = 300
    smooth = True
    smoothval = 3

    smooth1 = True
    smoothval1 = 2

    linearsmooth = False
    linearsmoothval = 0
    der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....
    der1 = ""

    poweraxis = ""
    zeropoweroffset = True
    calcpower = False
    zerocurve = True

    #trims data if nonzero
    xlimits = [-.5, 2]
    ylimits = [0]
    ylimits1 = [0]
    zlimit = [0,0] 
    # zlimit = [0,3000]
    # zlimit = [-.1,3.1]  

    xlog = False
    ylog = False
    xabs = False    # Take absolute value of the data in x or y
    yabs = False

    ############################################# DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    if direrction == "y":
        gatecuts = np.linspace(np.min(xval), np.max(xval), slices)
    if direrction == "x":
        gatecuts = np.linspace(np.min(yval), np.max(yval), slices)

    if poweraxis == "x":
        if calcpower:
            pval = []
            for i in range(xlen):
                pval.append(np.mean(power[:,i]))
            pval = np.asfarray(pval)
            xval = pval
        else:
            xval = xval - xval[0]
        if zeropoweroffset:
            offset_ = np.mean(data[:,0])
            data = data - offset_
    if poweraxis == "y":
        if calcpower:
            pval = []
            for i in range(ylen):
                pval.append(np.mean(power[i,:]))
            pval = np.asfarray(pval)
            yval = pval - pval[0]
        else:
            yval = yval - yval[0]
        if zeropoweroffset:
            offset_ = np.mean(data[0,:])
            data = data - offset_



    if np.min(xval) == np.max(xval):
        xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
    if np.min(yval) == np.max(yval):
        yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

    if smooth:
        data = applyGauss(data, smoothval)

    if linearsmooth:
        for i in range(ylen):
            data[i,:] = applyGauss(data[i,:], linearsmoothval)

    if interpolate:
        f = interp.interp2d(xval, yval, data, kind = 'linear')
        xval = np.linspace(np.min(xval), np.max(xval), newbox)
        yval = np.linspace(np.min(yval), np.max(yval), newbox)
        data = f(xval, yval)

    xlen = xval.size
    ylen = yval.size

    if der == 'x':
        derstep = np.abs(xval[0] - xval[1])
        ydelta, data = np.gradient(data, derstep)
        data = data * -1

    # if xlog:
    #     data = np.abs(data)


    ############################################# DATA END STUFF

    #Loads and builds axes
    axs, figs = customPlot(4.83, 2, sidebar = False)
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))


    # ax[0].yaxis.label.set_color('red')
    # ax[0].tick_params(axis='y', colors='red')

    ############################################# PLOT STUFF


    # mpl.colorbar.ColorbarBase(ax[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')

    data = data + offset

    if len(ax) > 1:
        ax[1].yaxis.tick_right()
        ax[1].yaxis.set_label_position('right')
        gatenorm = mpl.colors.Normalize(vmin = np.max(gatecuts), vmax = np.min(gatecuts))
        mpl.colorbar.ColorbarBase(ax[1], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='vertical')

    x = []
    y = []

    for i in range(len(gatecuts)):
        if direrction == "x":
            gateindex = np.searchsorted(yval, gatecuts[i])
            x.append(gatecuts[i])
            y.append(np.max((data[gateindex, :])))

        else:
            gateindex = np.searchsorted(xval, gatecuts[i])
            x.append(gatecuts[i])
            y.append(np.max((data[ :, gateindex])))


    for i in range(len(y)):
        y[i] = y[i]/72400

    ax[0].scatter(x, y, s = width, c = clr, alpha = alpha)


    # pl.figure()
    # pl.imshow(data)
    # pl.show()


    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    for i in vlines:
        ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits == [0]:
        if direrction == "x":
            ax[0].set_xlim(np.min(yval), np.max(yval))
        else:
            ax[0].set_xlim(np.min(xval), np.max(xval))
    else:
        ax[0].set_xlim(xlimits[0], xlimits[1])

    if ylimits != [0]:
        ax[0].set_ylim(ylimits[0], ylimits[1])


    #Sets log
    if xlog:
        ax[0].set_yscale("log")

    today = date.today()
    savename = name + "_customscatter_" + direrction +  str(today)

    # xt = np.linspace(0, .15, 6)
    # xt = np.around(xt, 2)

    # ax[0].set_xticks(xt, xt)

    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)


    if save:
        pl.savefig("C:/QMO/generated plots/Figs/" + savename + "_" + str(savemod) +".png")

    if show:
        pl.show()
    
    pl.clf()


# plines = []
# plines.append([1,0,0,0,0,0,0,0,0,0])
# plines.append([1,1,0,0,0,0,0,0,0,0])
# plines.append([1,1,1,0,0,0,0,0,0,0])
# plines.append([1,1,1,1,0,0,0,0,0,0])
# plines.append([1,1,1,1,1,0,0,0,0,0])
# plines.append([1,1,1,1,1,1,0,0,0,0])
# plines.append([1,1,1,1,1,1,1,0,0,0])
# plines.append([1,1,1,1,1,1,1,1,0,0])
# plines.append([1,1,1,1,1,1,1,1,1,0])
# plines.append([1,1,1,1,1,1,1,1,1,1])

# for i in range(len(plines)):
#     main(10+i, plines[i])

if __name__ == '__main__':
    main()