'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

# import matplotlib as mpl
# import matplotlib.pyplot as pl
# import matplotlib.colors as colors
# from matplotlib import cm
from msilib.schema import Directory
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

# name = "2022_02_03_12_alphamap" #Alpha fit maps
# name = "2022_02_03_12_errormap"

# name = "CPS_2022_02_01_17"  #Hi def bright current, wide range
# name = "2022_02_03_15_powercube"  #high power sd-gate powercube
# name1 = "CPS_2022_02_01_6"  #dark current
def main(savemod = 0, imagemod = 14):
    name = "2022_01_29_powermap"  #power, SD, hypercube at -5.13V gate
    name = "2022_01_29_alpha"

    path = "C:/QMO/Built"

    # xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
    # xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
    # xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)

    # xval, yval, cval, zval, data, rf = loader.load_hypercube(path, name)
    xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)

    # data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
    # data = data[:,:,33]

    # data = data - np.mean(data[:,:,0,14])
    # data = data[:,:,cval.size-2, imagemod]
    # rf = rf[:,:,cval.size-2, imagemod]

    data = data[:,:,1]
    cutoff = 8e6
    savemod = 103

    xrange = np.abs(xval[-1] - xval[0])
    xval = np.linspace(-xrange/2, xrange/2, xval.size)
    yrange = np.abs(yval[-1] - yval[0])
    yval = np.linspace(-yrange/2, yrange/2, yval.size)

    #4.8
    # yoff = 0
    # yoff = -4.5
    yoff = 0
    # xlimits = [-14, 14]
    # ylimits = [-3.39+yoff, 14 + yoff]
    xlimits = [-12, 9]
    ylimits = [-9, 11]
    xlabel = "$\it{x}$  (micron)"
    ylabel = "$\it{y}$  (micron)"
    xlabel = ""
    ylabel = ""
    # xlabel = ""
    zlabel = "$\it{I}$  (nA)"
    zlabel = ""

    # data = data * 1e9       # calibrates to real values

    save = True
    show = True

    #Parameters
    mapcolor_name = 'magma' #Color of map
    zerocolor = False
    contours = False

    hlines = []
    vlines = []

    #simple parameters -- in the otf version these pretty much stay the way they are
    interpolate = False
    newbox = 300
    smooth = True
    smoothval = 0

    linearsmooth = False
    linearsmoothval = 0
    der = ""        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

    poweraxis = ""

    #trims data if nonzero

    zlimit = [0,3000] 

    zlog = False
    ylog = False
    xabs = False    # Take absolute value of the data in x or y
    yabs = False

    ############################################# DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    
    if cutoff > 0:
        for i in range(xlen):
            for p in range(ylen):
                if data[p,i] > cutoff or data[p,i] < 0:
                    data[p,i] = 0


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
    if der == 'y':
        derstep = np.abs(yval[0] - yval[1])
        data, xdata = np.gradient(data, derstep)
        data = data * -1

    if zlog:
        data = np.abs(data)

    ############################################# DATA END STUFF

    #Loads and builds axes
    axs, figs = customPlot(4.83, 3, sidebar = True)
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))


    ############################################# PLOT STUFF

    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = 1)
    for i in vlines:
        ax[0].axvline(i, linestyle = ":", c = 'k', linewidth = 1)

    if zerocolor:
        zmax, zmin = np.max(data), np.min(data)
        if np.abs(zmax) > np.abs(zmin):
            zmin = -zmax
        else:
            zmax = -zmin

    cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
    xt = get_extent(xval, yval, inverty= True)

    ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position('right')
    ax[1].set_ylabel(zlabel)
    mpl.colorbar.ColorbarBase(ax[1], cmap = cmap, norm = cnorm, orientation='vertical')

    if contours:
        ax[0].contour( rf, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')

    #Sets edges of image
    if xlimits == [0]:
        None
    else:
        ax[0].set_xlim(xlimits[0], xlimits[1])

    if ylimits == [0]:
        None
    else:
        ax[0].set_ylim(ylimits[0], ylimits[1])


    today = date.today()
    savename = name + "_customimage_" + str(today)

    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)

    if save:
        pl.savefig("C:/QMO/generated plots/Figs/" + savename + "_" + str(savemod) +".png")

    if show:
        pl.show()

    pl.clf()

# for i in range(30):
#     main(i, i)

if __name__ == '__main__':
    main()