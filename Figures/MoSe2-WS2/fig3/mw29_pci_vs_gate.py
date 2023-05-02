'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib.font_manager as font_manager
import matplotlib as mpl
font_dir = ["C:/Helvetica World"]
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)
mpl.rcParams['font.family'] = 'Helvetica World'
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt, savgol_filter
# from curves import diode, exponential, twobody, ln, power, line

from datetime import date

import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader

def main(named = 0):

    name = "2021_03_25_51_powercube" #mw29 powermap
    name1 = "CPS_2021_03_25_29" #mw29 reverse bias

    path = "C:/QMO/Built"
    savename = "MW29_pci_v_gate"

    if named != 0:
        name = named

    # xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
    xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)

    # data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
    data = data[:,:,20] - np.mean(data[:,0:10,0])
    data = data * 1e9       # calibrates to real values
    # data = data + 2.137

    xlabel = "$V_{Gate}$  (V)"
    # zlabel = 
    ylabel = "$\it{I}$  (nA)"

    save = True
    show = True

    # mapcolor_name = 'seismic' #Color of map
    mapcolor_name = 'magma' #Color of map
    mapcolor_name = "plasma"
    zerocolor = False
    contours = False

    invertyaxis = True
    upper = "lower"

    hlines = []
    vlines = []

    #simple parameters -- in the otf version these pretty much stay the way they are
    interpolate = False
    newbox = 400
    smooth = 1
    savgolsmooth = 0
    savgolpoly = 3

    linearsmooth = False
    linearsmoothval = 0
    der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

    poweraxis = ""

    #trims data if nonzero

    xlimits = [-8,-2]
    ylimits = []
    # xlimits = [0]
    # ylimits = [0]
    zlimit = [1.6, 2, 2.4]
    zlimit = [2.4, 2, 1.6]
    clr = ['k', 'b', 'r']
    clr = []
    colorpadding = 2

    zlog = False
    ylog = False
    xabs = False    # Take absolute value of the data in x or y
    yabs = False

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    # for i in range(xlen):
    #     for p in range(ylen):
    #         if data[p,i] == np.inf:
    #             data[p,i] = 0


    if np.min(xval) == np.max(xval):
        xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
    if np.min(yval) == np.max(yval):
        yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

    if smooth > 0:
        data = applyGauss(data, smooth)

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
    if der == 'y':
        derstep = np.abs(yval[0] - yval[1])
        data, xdata = np.gradient(data, derstep)
        data = data * -1

    if zlog:
        # data = np.sign(data) * np.log(np.abs(data + 1))
        # data = np.sign(data)
        data = np.sign(data) * np.log(np.abs(data)+1)
        # data = np.log(data)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA END STUFF

    #Loads and builds axes
    axs, figs = customPlot(4, 4, sidebar = False)
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT STUFF

    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .8)
    for i in vlines:
        ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 2, alpha = .5)


    if clr == []:
        clr = color_meline(mapcolor_name, len(zlimit) + colorpadding)

    x = yval
    for i in range(len(zlimit)):
        zindex = np.searchsorted(xval, zlimit[i])
        y = data[:,zindex]
        if savgolsmooth > 1:
            y = savgol_filter(y, savgolsmooth, savgolpoly)

        ax[0].plot(x, y, linewidth = 2, c = clr[i])



    # ax[1].yaxis.tick_right()
    # ax[1].yaxis.set_label_position('right')
    # mpl.colorbar.ColorbarBase(ax[1], cmap = cmap, norm = cnorm, orientation='vertical')
    # ax[1].set_ylabel(zlabel)

    if contours:
        ax[0].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

    #Sets edges of image
    if xlimits == []:
        None
    else:
        ax[0].set_xlim(xlimits[0], xlimits[1])

    if ylimits == []:
        None
    else:
        ax[0].set_ylim(ylimits[0], ylimits[1])


    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save

    if save:
        pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

    if show:
        pl.show()

if __name__ == '__main__':
    main()