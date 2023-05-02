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
from curves import simplepower

def main(savemod = 0, plotgline = 0):

    # name = "CPS_2021_03_18_1"
    # name = "CPS_2021_11_15_18_filt"
    # name = "CPS_2021_11_12_42"
    # name = "CPS_2021_12_09_37"
    # name1 = "CPS_2022_02_01_6"
    # name = "2022_02_03_12_alphamap"
    # name = "2022_02_03_12_errormap"
    name = "2021_03_25_51_powercube"
    name1 = "CPS_2022_02_01_6"
    # name = "CPS_2022_01_26_55"
    # name = "2022_02_03_12_powercube"
    # name1 = "2022_02_03_15_powercube"

    double = False

    # name = "CPS_2022_02_01_6"

    path = "C:/QMO/Built"

    # xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
    # xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
    xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
    if double:
        # xval1, yval1, cval1, data1, rf1, power1, info1 = loader.load_cube(path, name1)
        # xval1, yval1, data1, rf1, info1 = loader.load_map(path, name1, returnpower=False)
        xval1, yval1, data1, rf1, power1, info1 = loader.load_map(path, name1, returnpower=True)

    # data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
    # data = data[:,:,33]


    xlabel = ""#"$V_{TG} = V_{BG}$  (V)"
    ylabel = ""#"$\it{I}$  (nA)"
    ylabel1 = ""#"$\it{I}$  (nA)"

    data = data * 1e9       # calibrates to real values
    if double:
        data1 = data1 * 1e9

    # data = np.flip(data, 1)
    # data = data * -1
    # if double:
    #     data1 = data1 * -1
    # xval = np.flip(xval)
    # data = np.flip(data, 1)

    # data1 = np.flip(data1, 0)
    # yval1 = np.flip(yval1)
    # xval1 = np.flip(xval1)
    # data1 = np.flip(data1, 1)

    save = True
    show = True

    savemod = savemod

    #Parameters
    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices
    contours = False

    colormap_name = "viridis_r"

    #How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
    slices = 50
    # gatecuts = np.linspace(-5, -6, 10)
    setgate = -5.2
    gatecuts = np.linspace(0, 9, 7)
    gatecutsindex = [0, 3, 6, 9, 12, 15, 18]
    # gatecuts = np.asarray([-1, -2, -3, -3.5, -4, -4.5, -5])
    # gatecuts = [-.5]
    # plotgline = [1, 1, 1, 1, 1, 1, 1]
    plotgline = np.ones(gatecuts.shape)
    # plotgline = plotgline
    speccolor = color_meline(colormap_name, gatecuts.size  )
    # speccolor = ["r", 'b', 'g']
    style = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
    alpha = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    width = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

    gatecuts1 = [-5.2]
    plotgline1 = [1]
    speccolor1 = ["k"]
    style1 = ["-"]
    alpha1 = [1]
    width1 = [2]

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
    smoothval = 1

    smooth1 = True
    smoothval1 = 1

    linearsmooth = False
    linearsmoothval = 0
    der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....
    der1 = ""

    poweraxis = "y"
    zeropoweroffset = True

    #trims data if nonzero
    xlimits = [.5, 2]
    ylimits = [-1.6, .5]
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
        yval = pval - pval[0]
        if zeropoweroffset:
            offset_ = np.mean(data[0,:])
            data = data - offset_

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
    axs, figs = customPlot(4.83, 3, sidebar = True)
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))


    # ax[0].yaxis.label.set_color('red')
    # ax[0].tick_params(axis='y', colors='red')

    if double:
        ax1 = ax[0].twinx()

    ############################################# PLOT STUFF

    if gatecuts != []:
        if speccolor == []:
            clrgate = color_meline(colormap_name, len(gatecuts)  )
        else:
            clrgate = speccolor


    # mpl.colorbar.ColorbarBase(ax[0], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='horizontal')

    if len(ax) > 1:
        ax[1].yaxis.tick_right()
        ax[1].yaxis.set_label_position('right')
        gatenorm = mpl.colors.Normalize(vmin = np.max(gatecuts), vmax = np.min(gatecuts))
        mpl.colorbar.ColorbarBase(ax[1], cmap = pl.get_cmap(colormap_name), norm = gatenorm, orientation='vertical')

    for i in range(len(gatecuts)):
        if direrction == "x":
            if gatecutsindex == []:
                gateindex = np.searchsorted(yval, gatecuts[i])
            else:
                gateindex = gatecutsindex[i]
            print(gateindex)
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

    if double:
        for i in range(len(gatecuts1)):
            if plotgline1[i] == 1:
                if direrction == "x":
                    gateindex = np.searchsorted(yval1, gatecuts1[i])
                    x = xval1
                    y = (data1[gateindex, :] + offset1) * scaling1
                else:
                    gateindex = np.searchsorted(xval1, gatecuts1[i])
                    x = yval1
                    y = (data1[:, gateindex] + offset1) * scaling1

                if plotgline1[i] == 1:
                    ax1.plot(x, y, linewidth = width1[i], c = speccolor1[i], linestyle = style1[i], alpha = alpha1[i])


    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    for i in vlines:
        ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits == [0]:
        if direrction == "x":
            ax[0].set_xlim(np.min(xval), np.max(xval))
        else:
            ax[0].set_xlim(np.min(yval), np.max(yval))
    else:
        ax[0].set_xlim(xlimits[0], xlimits[1])

    if ylimits != [0]:
        ax[0].set_ylim(ylimits[0], ylimits[1])
    if double:
        if ylimits1 != [0]:
            ax1.set_ylim(ylimits1[0], ylimits1[1])

    #Sets log
    if xlog:
        ax[0].set_yscale("log")

    today = date.today()
    savename = name + "_custom_" + direrction +  str(today)

    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    if double:
        ax1.set_ylabel(ylabel1)

    if save:
        pl.savefig("C:/QMO/generated plots/Figs/" + savename + "_" + str(savemod) +".png")

    if show:
        pl.show()
    
    pl.clf()


plines = []
plines.append([1,0,0,0,0,0,0,0,0,0])
plines.append([1,1,0,0,0,0,0,0,0,0])
plines.append([1,1,1,0,0,0,0,0,0,0])
plines.append([1,1,1,1,0,0,0,0,0,0])
plines.append([1,1,1,1,1,0,0,0,0,0])
plines.append([1,1,1,1,1,1,0,0,0,0])
plines.append([1,1,1,1,1,1,1,0,0,0])
# plines.append([1,1,1,1,1,1,1,1,0,0])
# plines.append([1,1,1,1,1,1,1,1,1,0])
# plines.append([1,1,1,1,1,1,1,1,1,1])

# for i in range(len(plines)):
#     main(10+i, plines[i])

if __name__ == '__main__':
    main()