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
from curves import simplepower, twobody1, twobody

def main(savemod = 0, plotgline = 0):

    name = "CPS_2022_01_26_53"
    name1 = "CPS_2022_01_26_55"

    double = True
    fitfunc = twobody
    trimval = 6

    # name = "CPS_2022_02_01_6"

    path = "C:/QMO/Built"

 
    xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
    if double:
        xval1, yval1, data1, rf1, power1, info1 = loader.load_map(path, name1, returnpower=True)

    xlabel = ""#"$V_{TG} = V_{BG}$  (V)"
    ylabel = ""#"$\it{I}$  (nA)"
    ylabel1 = ""#"$\it{I}$  (nA)"

    data = data * 1e9       # calibrates to real values
    data1 = data1 * 1e9       # calibrates to real values

    save = True
    show = True

    savemod = 0

    #Parameters
    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices
    contours = False

    colormap_name = "viridis_r"

    #How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
    slices = 50

    gatecuts = [10]
    gatecutsindex = [0, 3, 6, 9, 12, 15, 18]
    plotgline = [1, 1, 1, 1, 1, 1, 1]

    speccolor = ["k", 'b', 'k']
    speccolor = []
    style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
    alpha = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    width = [ 2, 1.6, 2, 2, 2, 2, 2, 2, 2]

    gatecuts1 = [ -5.13]
    plotgline1 = [1, 0, 0, 1]
    speccolor1 = ["k", 'b', 'k']
    style1 = ["-", "-", "-", "-"]
    alpha1 = [1, 1, 1, 1]
    width1 = [ 2, 1.6, 2]

    scaling = 1
    scaling1 = 1
    offset = 0
    offset1 = 0

    direrction = "x"

    hlines = [0]
    vlines = [0]

    #simple parameters -- in the otf version these pretty much stay the way they are
    interpolate = False
    newbox = 300
    smooth = True
    smoothval = 0

    smooth1 = True
    smoothval1 = 0

    linearsmooth = False
    linearsmoothval = 0
    der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....
    der1 = ""

    poweraxis = "y"
    zeropoweroffset = True

    #trims data if nonzero
    xlimits = [-.5, 3]
    ylimits = [0]
    ylimits1 = [-2, 40]
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
    if double:
        xlen1 = xval1.size
        ylen1 = yval1.size

    # for i in range(xlen):
    #     for p in range(ylen):
    #         if data[p,i] == np.inf:
    #             data[p,i] = 0

    if poweraxis == "x":
        pval = []
        for i in range(xlen):
            pval.append(np.mean(power[:,i]))
        pval = np.asfarray(pval)
        xval = pval
        if zeropoweroffset:
            offset_ = np.mean(data[:,0])
            data = data - offset_
    if poweraxis == "y":
        pval = []
        for i in range(ylen):
            pval.append(np.mean(power[i,:]))
        pval = np.asfarray(pval)
        yval = pval
        if zeropoweroffset:
            offset_ = np.mean(data[0,:])
            data = data - offset_
        
        pval1 = []
        for i in range(ylen1):
            pval1.append(np.mean(power1[i,:]))
        pval1 = np.asfarray(pval1)
        yval1 = pval1
        if zeropoweroffset:
            offset1_ = np.mean(data1[0,:])
            data1 = data1 - offset1_
        

    # pl.figure()
    # pl.imshow(data)
    # pl.show()

    # data = data + 21

    if np.min(xval) == np.max(xval):
        xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
    if np.min(yval) == np.max(yval):
        yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

    if smooth:
        data = applyGauss(data, smoothval)
    if smooth1 and double:
        data1 = applyGauss(data1, smoothval1)

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
    axs, figs = customBottomtop(5, 5, sidebar = False)
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))

    ax[0].set_xticks([])
    # ax[0].yaxis.label.set_color('red')
    # ax[0].tick_params(axis='y', colors='red')

    if double:
        ax1 = ax[1]

    ############################################# PLOT STUFF

    usegateindex = False
    if len(gatecutsindex) > 0:
        usegateindex = True
        gatecuts = gatecutsindex

    if gatecuts != []:
        if speccolor == []:
            clrgate = color_meline(colormap_name, len(gatecuts)  )
        else:
            clrgate = speccolor


    # mpl.colorbar.ColorbarBase(ax[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')

    if len(ax) > 2:
        ax[2].yaxis.tick_right()
        ax[2].yaxis.set_label_position('right')
        gatenorm = mpl.colors.Normalize(vmin = np.max(gatecuts), vmax = np.min(gatecuts))
        mpl.colorbar.ColorbarBase(ax[2], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='vertical')



    for i in range(len(gatecuts)):
        if direrction == "x":
            if usegateindex:
                gateindex = gatecuts[i]
            else:
                gateindex = np.searchsorted(yval, gatecuts[i])
            x = xval
            y = (data[gateindex, :] + offset) * scaling
        else:
            gateindex = np.searchsorted(xval, gatecuts[i])
            x = yval
            y = (data[:, gateindex] + offset) * scaling

        if plotgline[i] == 1:
            ax[0].plot(x, y, linewidth = width[i], c = clrgate[i], linestyle = style[i], alpha = alpha[i])
        else:
            ax[0].plot(x, y, linewidth = width[i], c = clrgate[i], linestyle = style[i], alpha = 0)

    # pl.figure()
    # pl.imshow(data)
    # pl.show()

    ptrim = np.searchsorted(yval1, trimval)

    if double:
        alpha = []
        sdx = []
        for i in range(xlen1):
            x = yval1
            y = data1[:,i]

            try:
                par, pcov = curve_fit(fitfunc, x[:ptrim], y[:ptrim], p0 = [10,1], maxfev = 3200)
                error = np.sum( (np.abs(y[:ptrim]) - np.abs(fitfunc(x[:ptrim], *par))) )
                e0 = np.abs(pcov[0,0])
                # if error < 1.37:
                if e0 < 10:
                # if par[0] < 100 and par[0] > 0:
                    sdx.append(xval1[i])
                    alpha.append(par[0])
                else:
                    sdx.append(xval1[i])
                    alpha.append(0)


                # pl.figure()
                # pl.scatter(x[:ptrim], y[:ptrim], c = "k")
                # pl.plot(x[:ptrim], fitfunc(x[:ptrim], *par), c = "r")

            except RuntimeError:
                sdx.append(xval1[i])
                alpha.append(0)

        ax1.plot(sdx, alpha, linewidth = width1[0], c = speccolor1[0], linestyle = style1[0], alpha = alpha1[0])


    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
        ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    for i in vlines:
        ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 1)
        ax[1].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits == [0]:
        if direrction == "x":
            ax[0].set_xlim(np.min(xval), np.max(xval))
            ax[1].set_xlim(np.min(xval), np.max(xval))
        else:
            ax[0].set_xlim(np.min(yval), np.max(yval))
            ax[1].set_xlim(np.min(yval), np.max(yval))
    else:
        ax[0].set_xlim(xlimits[0], xlimits[1])
        ax[1].set_xlim(xlimits[0], xlimits[1])

    if ylimits != [0]:
        ax[0].set_ylim(ylimits[0], ylimits[1])
    if ylimits1 != [0]:
        ax1.set_ylim(ylimits1[0], ylimits1[1])

    #Sets log
    if xlog:
        ax[0].set_yscale("log")

    today = date.today()
    savename = name + "_custom_spec_" + direrction +  str(today)

    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax1.set_ylabel(ylabel1)

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