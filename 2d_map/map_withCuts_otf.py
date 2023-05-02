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

from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map

def main():
    year = 2022
    month = 1
    day = 26
    run = 55
    
    newpath = "C:/Jed/Built"
    path = "E:/Data/Raw/%i/%i_%.2i/%i_%.2i_%.2i/" % (year, year, month, year, month, day)
    name = 'CPS_%i_%.2i_%.2i_%i' % (year, month, day, run)
    build_map(path, name, newpath)

    path = "C:/Jed/Built"
    # xval, yval, data, rf, power, info = loader.load_map(path, name, returnpower=True)
    xval, yval, data, rf, info = loader.load_map(path, name, returnpower=False)
    info = make_dict(info)


    save = False
    savename = name + "_otf"
    show = True

    #Parameters
    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices

    #How many linecuts to take, evenly spaced. Can also take a list of gate/sd values and an on/off int switch
    slices = 50
    gateCuts = np.linspace(np.min(yval), np.max(yval), slices)
    # gateCuts = np.linspace(-3.5, -1, slices)
    plotgline = np.ones(slices)
    sourceCuts = np.linspace(np.min(xval), np.max(xval), slices)
    plotsdline = np.ones(slices)

    #simple parameters -- in the otf version these pretty much stay the way they are
    interpolate = False
    newbox = 1000
    smooth = True
    smoothval = 0
    der = ''        #take derivatives of data slices? 0 = no derivative, 1 = 1st derivative, 2 etc....

    xlimits = [0]
    x1limits = [0]
    ylimits = [0]
    zlimit = [0,0] 
    zlog = False
    invert_y = False
    xlog = False
    ylog = False
    xstep = 1   #-- how many line cutes to keep in each direction, "1" keeps all lines
    ystep = 1
    xabs = False    # Take absolute value of the data in x or y
    yabs = False
    makenegative = False #swtiches the sign of the data

    data = data * 1e9       # calibrates to real values
    # data = data + 0.20478921940480618

    ############################################# DATA DO STUFF

    if yval[1] < yval[0]:
        yval = np.flip(yval)
        data = np.flip(data, 0)
    if xval[1] < xval[0]:
        xval = np.flip(xval)
        data = np.flip(data, 1)

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
        data = np.abs(data)
    if yabs:
        ydata = np.abs(ydata)

    ############################################# DATA END STUFF

    #Loads and builds axes
    axs, axcbs, figs = makeaxes(3, cbar = True)
    ax = []
    fig = pl.figure(figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))
    axcb = fig.add_axes(axcbs[0])
    axcb.set_xticks([])
    axcb.yaxis.tick_right()
    axcb.yaxis.set_label_position('right')

    ############################################# PLOT STUFF

    no = len(gateCuts)
    if gateCuts != []:
        clr = color_meline(colormap_name, no)    

    if zlimit == [0]:
        zlimit = [0 , 0]
    
    # if zlog:
    #     cmap, cnorm = color_memap(mapcolor_name, np.log(np.abs(data)) , dmin = zlimit[0], dmax = zlimit[1])
    # else:
    cmap, cnorm = color_memap(mapcolor_name, data , dmin = zlimit[0], dmax = zlimit[1])

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

    # if zlog:
    #     ax[0].imshow(np.log(np.abs(data)) , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
    # else:
    ax[0].imshow(data , cmap = cmap, norm = cnorm, extent = xt, aspect = 'auto', origin = 'lower')
    mpl.colorbar.ColorbarBase(axcb, cmap = cmap, norm = cnorm, orientation='vertical')

    if gateCuts != []:
        clrgate = color_meline(colormap_name, len(gateCuts) + 3 )
    if sourceCuts != []:
        clrsd = color_meline(colormap_name, len(sourceCuts) + 3 )

    ymin1 = 10
    ymax1 = -10
    ymin2 = 10
    ymax2 = -10

    for i in range(len(gateCuts)):
        if plotgline[i] == 1:
            None
            # ax[0].axhline(gateCuts[i], linewidth = 1, linestyle = ':', c = clrgate[i])
        gateindex = np.searchsorted(yval, gateCuts[i])
        x = xval
        y = data[gateindex, :]
        if plotgline[i] == 1:
            if zlog:
                y = np.abs(y)
            ax[1].plot(x, y, linewidth = 1.5, c = clrgate[i])

        #Gets min and maxes of data so that the script can scale properly
        if np.min(y) < ymin1:
            ymin1 = np.min(y)
        if np.max(y) > ymax1:
            ymax1 = np.max(y)


    for i in range(len(sourceCuts)):
        if plotsdline[i] == 1:
            None
            # ax[0].axhline(gateCuts[i], linewidth = 1, linestyle = ':', c = clrgate[i])
        sdindex = np.searchsorted(xval, sourceCuts[i])
        x = yval
        y = data[: , sdindex]
        if plotsdline[i] == 1:
            if zlog:
                y = np.abs(y)
            ax[2].plot(x, y, linewidth = 1.5, c = clrsd[i])

        #Gets min and maxes of data so that the script can scale properly
        if np.min(y) < ymin2:
            ymin2 = np.min(y)
        if np.max(y) > ymax2:
            ymax2 = np.max(y)
        

    ymaxrange1 = ymax1 - ymin1
    ymaxrange2 = ymax2 - ymin2
    margin1 = np.abs(ymaxrange1) * .05
    margin2 = np.abs(ymaxrange2) * .05
    
    if invert_y:
        ax[1].set_ylim( ymax1 + margin1 , ymin1 - margin1 )
        ax[2].set_ylim( ymax2 + margin2 , ymin2 - margin2 )
    else:
        ax[1].set_ylim( ymin1 - margin1 , ymax1 + margin1 )
        ax[2].set_ylim( ymin2 - margin2 , ymax2 + margin2 )

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits == [0]:
        ax[1].set_xlim(np.min(xval), np.max(xval))
        ax[2].set_xlim(np.min(yval), np.max(yval))
        # ax[1].set_xticks( np.linspace(np.min(xval), np.max(xval), 5) )
    else:
        ax[1].set_xlim(x1limits[0], x1limits[1])
        ax[1].set_xticks( xtick )
        ax[0].set_xticks( xtick )

    if ylimits != [0]:
        ax[1].set_ylim(ylimits[0], ylimits[1])

    # ax[1].set_yticks( yt )
    # ax[0].set_yticks( ytick )
    # ax[1].set_yticks( ytick1 )

    #Sets log
    if zlog:
        ax[1].set_yscale('log')
    if xlog:
        ax[0].set_yscale("log")
    if ylog:
        ax[1].set_yscale("log")
    

    ############################################# END PLOT STUFF
    ## Labels

    fig.suptitle(info["Run Number"])
    ax[0].set_title("Map")
    ax[1].set_title("Horizontal Cuts")
    ax[2].set_title("Vertical Cuts")

    ax[0].set_xlabel(info["Fast Axis Variable"])
    ax[0].set_ylabel(info["Slow Axis Variable"])

    ax[1].set_xlabel(info["Fast Axis Variable"])
    ax[1].set_ylabel("Current (nA)")

    ax[2].set_xlabel(info["Slow Axis Variable"])
    ax[2].set_ylabel("Current (nA)")

    if interpolate:
        interpbool = "interpolated"
    else:
        interpbool = ""

    fig.suptitle(name + "  :  gaussianFilter = " + str(smoothval) + "  :  " + interpbool, fontsize = 10)

    zlabel = "Current (nA)"
    if "der" == "x":
        zlabel = "Conduction (nA / V)"

    axcb.set_ylabel(zlabel)

    if save:
        pl.savefig("C:/Jed/Figures/" + savename + ".png")
    
    if show:
        pl.show()


if __name__ == '__main__':
    main()
