'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from scipy import ndimage
from scipy import interpolate as interp

import sys
sys.path.append("C:/Users/Jedki/Documents/analCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map

def main():
    oldpath = "C:/Users/Jedki/QMOdata/Raw Data/2021_10/2021_10_19"
    newpath = "C:/Users/Jedki/QMOdata/Built"
    name = "CPS_2021_10_19_6"
    build_map(oldpath, name, newpath)

    path = "C:/Users/Jedki/QMOdata/Built"
    xval, yval, data, rf, info = loader.load_map(path, name)
    info = make_dict(info)

    # yval = np.flip(yval)
    # data = np.flip(data, 0)

    #Plot information
    title1 = "photocurrent"
    title2 = "reflection"
    x1label = 'x'
    y1label = 'y'
    zlabel = 'Current (nA)'
    x2label = 'x'

    y2label = ""


    save = True
    show = True
    contour = True
    crosshairs = True

    savename = name + "_otf"
    if contour:
        savename = savename + "_rfcont"

    #Parameters
    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'gray_r' #Color of slices

    logfile = oldpath + "/" + name + "_log.log"
    f = open(logfile, "r")
    info = {}
    info["scan"] = name

    for line in f:
        s = line.split(":")
        info[s[0]] = s[1]
    f.close()

    interpolate = False
    newbox = 1000
    smooth = True
    smoothval = 0
    der = ''        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

    xlimits = [0]
    x1limits = [0]
    ylimits = [0]
    zlimit = [0, 0] 
    invert_y = False
    xstep = 1   #-- how many line cutes to keep in each direction, "1" keeps all lines
    ystep = 1
    xabs = False    # Take absolute value of the data in x or y
    yabs = False
    makenegative = False #swtiches the sign of the data

    data = data * 1e9       # calibrates to real values

    ############################################# DATA DO STUFF

    if makenegative:
        data = data * -1

    if der == 'x':
        ydelta, data = np.gradient(data)
    elif der == 'y': 
        data, xdelta = np.gradient(data)

    if interpolate:
        f = interp.interp2d(xval, yval, data, kind = 'linear')
        xval = np.linspace(np.min(xval), np.max(xval), newbox)
        yval = np.linspace(np.min(yval), np.max(yval), newbox)
        data = f(xval, yval)

    if smooth:
        data = applyGauss(data, smoothval)

    # Makes data absolute
    if xabs:
        xdata = np.abs(xdata)
    if yabs:
        ydata = np.abs(ydata)

    ############################################# DATA END STUFF

    #Loads and builds axes
    axs, axcbs, figs = makeaxes(2, cbar = True)
    ax = []
    fig = pl.figure(figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))
    axcb = fig.add_axes(axcbs[0])
    axcb.set_xticks([])
    axcb.yaxis.tick_right()
    axcb.yaxis.set_label_position('right')

    ############################################# PLOT STUFF
   
    if zlimit == [0]:
        zlimit = [0 , 0]
    
    cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])
    rcmap, rcnorm = color_memap(colormap_name, rf , dmin = 0, dmax = 0)

    if xlimits == [0]:
        xr = xval
    else:
        xr = xval[np.searchsorted(xval, xlimits[0]) : np.searchsorted(xval, xlimits[1]) ]
        data = data[ : , np.searchsorted(xval, xlimits[0]) : np.searchsorted(xval, xlimits[1]) ]
        xval = xval [np.searchsorted(xval, xlimits[0]) : np.searchsorted(xval, xlimits[1]) ]

    if ylimits == [0]:
        yr = yval
    else:
        yr = yval[np.searchsorted(yval, ylimits[0]) : np.searchsorted(yval, ylimits[1]) ]

    xt = get_extent(xr, yr, inverty= True)

    ax[1].imshow(rf , cmap = rcmap, norm = rcnorm, extent = xt, aspect = 'auto')
    ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto')
    if contour:
        ax[0].contour(rf , cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'upper')
    if crosshairs:
        ax[1].axhline(np.mean(yr), c = 'r')
        ax[0].axhline(np.mean(yr), c = 'r')
        ax[1].axvline(np.mean(xr), c = 'r')
        ax[0].axvline(np.mean(xr), c = 'r')
    mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')


    # ax[1].set_yticks( yt )
    # ax[0].set_yticks( ytick )
    # ax[1].set_yticks( ytick1 )


    ############################################# END PLOT STUFF
    ## Labels

    ax[0].set_title(title1)
    ax[1].set_title(title2)


    ax[0].set_xlabel(x1label)
    ax[0].set_ylabel(y1label)

    ax[1].set_xlabel(x2label)
    ax[1].set_ylabel(y2label)


    if interpolate:
        interpbool = "interpolated"
    else:
        interpbool = ""

    # fig.suptitle(name + "  ||   "+ str(float(info["Filtered Wavelength"])) + "nm    ||   Gate = " + str(float(info["Backgate Start"])) + "V  Source = " + str(float(info["Source/Drain Start"])) +  "V   ||    gaussianFilter = " + str(smoothval) , fontsize = 10)
    fig.suptitle(name + "    ||   Gate = " + str(float(info["Backgate Start"])) + "V  Source = " + str(float(info["Source/Drain Start"])) +  "V   ||    gaussianFilter = " + str(smoothval) , fontsize = 10)

    axcb.set_ylabel(zlabel)

    if save:
        pl.savefig("C:/Users/jedki/QMOdata/generated plots/" + savename + ".png")
    
    if show:
        pl.show()


if __name__ == '__main__':
    main()
