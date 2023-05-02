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
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map
from curves import simplepower, twobody
from custom_maps import *

gates = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
gates = [14]

for gateindex in gates:

    axs, figs = fig_4()
    ax = []
    fig = pl.figure(num = 0, figsize=figs)
    for a in axs:
        ax.append(fig.add_axes(a))
    ax[3].set_xticks([])
    ax[3].set_yticks([])

    savename = "fig_4_full_built_" + str(gateindex)

    save = True
    show = False


    '''
                                                sd - power
    '''

    name = "2021_03_25_51_powercube"
    path = "C:/QMO/Built"


    xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
    data = data * 1e12       # calibrates to real values
    dataoffset = .0103+9.84
    # data = data - np.mean( data[yval.size - 20:yval.size+20,xval.size-20:xval.size+20,0] )
    # data = data[:,:,33]

    xlabel = "$V_{SD}$ (V)"
    xlabel = ""
    ylabel = "$\itI_{PC}$  (pA)"
    zlabel = r'Power $(\mu W)$'
    fontsize = 14
    ticksize = 14
    linewidth = 5
    markersize = 15



    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices

    gate = -6.5
    powerslice = [0, .01, .02, .03] # Go check cval after power arranging
    speccolor = []

    hlines = []
    vlines = [1.047]
    vlcolor = 'r'
    vlines = [1., 1.5]
    va = .85

    fill = []
    fill = []
    gradfill = [.5, 2.4]
    fillclr = ['deepskyblue', 'darkorange']
    fillclr = ['red']
    filla = .1

    smooth = 2
    savgolsmooth = 0
    savgolpoly = 3

    poweraxis = "z"
    zeropoweroffset = True

    #trims data if nonzero
    xlimits = [.25, 2.25]
    xticks = [.5, 1, 1.5, 2]
    ylimits = [-1,160]


    ############################################# DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    if poweraxis == "z":
        if zeropoweroffset:
            offset_ = np.mean(data[:,:,0])
            data = data - offset_ + dataoffset
            cval = cval - np.min(cval)

    if smooth > 0:
        for i in range(cval.size):
            data[:,:,i] = ndimage.gaussian_filter(data[:,:,i], smooth)

    powerslice = cval[:7:1]
    zticks = [0, 4, 8]

    ############################################# DATA END STUFF

    ############################################# PLOT STUFF


    # colorxt = [gradfill[0], gradfill[1], ylimits[0], ylimits[1]]
    # colordata = np.asarray([np.linspace(0, 1, 1000), np.linspace(0, 1, 1000)])
    # cmap, cnorm = custom_color(0, 1)
    # ax[0].imshow(colordata, cmap = cmap, norm = mpl.colors.Normalize(vmin = 0, vmax = 2), extent = colorxt, aspect = 'auto', origin = 'upper', alpha = .2)

    if speccolor == []:
        clr = color_meline(colormap_name, len(powerslice))
    else:
        clr = speccolor

    colordata = np.asarray([np.linspace(powerslice[0]*1000, powerslice[-1]*1000, 1000), np.linspace(powerslice[0]*1000, powerslice[-1]*1000, 1000)])
    colorxt=[powerslice[0]*1000, powerslice[-1]*1000, -1, 1]
    ax[-1].imshow(colordata, cmap = pl.get_cmap('viridis_r'), norm = mpl.colors.Normalize(vmin = powerslice[0]*1000, vmax = powerslice[-1]*1000), extent = colorxt, aspect = 'auto', origin = 'upper')
    ax[-1].set_yticks([])
    # ax[-1].yaxis.tick_right()
    # ax[-1].yaxis.set_label_position('right')
    ax[-1].set_xticks(zticks, color = 'w')

    x = xval
    for i in range(len(powerslice)):
        powerindex = np.searchsorted(cval, powerslice[i])
        # gateindex = np.searchsorted(yval, gate)
        y = data[gateindex, :, powerindex]

        if savgolsmooth > 1:
            y = savgol_filter(y, savgolsmooth, savgolpoly)

        # if smooth > 0:
        #     y= applyGauss(y, smooth)

        ax[0].plot(x, y, linewidth = linewidth, c = clr[i])


    for i in hlines:
        ax[0].axhline(i, linestyle = ":", c = 'k', linewidth = linewidth)
    for i in vlines:
        if i == 1:
            ax[0].axvline(i, ymin = 0, ymax = .5, linestyle = ":", c = vlcolor, linewidth = linewidth, alpha = va)
        else:
            ax[0].axvline(i, ymin = 0, ymax = .8, linestyle = ":", c = vlcolor, linewidth = linewidth, alpha = va)
        # ax[0].axvline(i, linestyle = "-", c = 'k', linewidth = 2, alpha = .25)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits != []:
        ax[0].set_xlim(xlimits[0], xlimits[1])
    if ylimits != []:
        ax[0].set_ylim(ylimits[0], ylimits[1])

    ax[0].set_xticks(xticks, ["", '', '', ''])

    # for i in range(len(fill)-1):
    #     ax[0].axvspan(fill[i], fill[i+1], color = fillclr[i], alpha = filla)

    ax[0].tick_params(axis='x', labelsize = ticksize, direction = 'out')
    ax[0].tick_params(axis='y', labelsize = ticksize, direction = 'in')
    ax[-1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'w')
    ax[0].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
    ax[0].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)
    ax[-1].set_xlabel(zlabel, fontsize = fontsize, labelpad = 0, rotation = 0)


    '''
                                                            Power scatter
    '''

    name = "2021_03_25_51_powercube"

    fit = True
    fitfunc = twobody
    trimval = 150

    maxval = True
    path = "C:/QMO/Built"


    xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)
    data = data * 1e12       # calibrates to real values
    cval = cval * 1000

    xlabel = r'Laser Power ($\mu $W)'
    ylabel = "$\itI_{PC}$ (pA)"
    fontsize = 14
    ticksize = 14
    linewidth = 3
    markersize = 20

    hlines1 = []
    vlines1 = []

    smooth = 1
    linearsmooth = 0
    savgolsmooth = 0
    savgolpoly = 3

    gate = [-6.5]
    source = [1.15]
    clr = ['k']
    lineclr = ['r']

    poweraxis = "z"
    zeropoweroffset = True
    calcpower = False
    zerocurve = True

    xlimits1 = [-.015, .2]
    ylimits = [-.03, .4]
    xlimits1 = [0, trimval]
    ylimits = [0, 120]
    xticks = [0, 50, 100]
    yticks = [.1, .2, .3, .4]
    yticks = [0, 50, 100]

    # xticks = [0, .5, ]

    ############################################# DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    if poweraxis == "z":
        if zeropoweroffset:
            offset_ = np.mean(data[:,:,0])
            data = data - offset_
            cval = cval - np.min(cval)


    if np.min(xval) == np.max(xval):
        xval = np.linspace(xval[0] - 1, xval[0] + 1, xlen)
    if np.min(yval) == np.max(yval):
        yval = np.linspace(yval[0] - 1, yval[0] + 1, ylen)

    ############################################# PLOT STUFF

    # data = data + offset

    x = cval
    for i in range(len(gate)):
        # gateindex = np.searchsorted(yval, gate[i])
        sourceindex = np.searchsorted(xval, source[i])

        y = data[gateindex, sourceindex, :]

        if zerocurve:
            y = y - y[0]

        if savgolsmooth > 1:
            y = savgol_filter(y, savgolsmooth, savgolpoly)

        if linearsmooth > 0:
            y= applyGauss(y, linearsmooth)

        ax[1].scatter(x, y, s = markersize, color = clr[i], alpha = 1, zorder = 2)

        if fit:
            try:
                print("Trying fit")
                ptrim = np.searchsorted(cval, trimval)
                par, pcov = curve_fit(fitfunc, x[:ptrim], y[:ptrim], maxfev = 3200)
                print("A = %f" % (par[0]))
                print("B = %f" % (par[1]))
                x_ = np.linspace(x[0], x[-1], 1000)
                ax[1].plot(x_, fitfunc(x_, *par), linewidth = linewidth, c = lineclr[i], zorder = 1)
            except RuntimeError:
                print("Fit failed")


    for i in hlines1:
        ax[1].axhline(i, linestyle = ":", c = 'k', linewidth = linewidth)
    for i in vlines1:
        ax[1].axvline(i, linestyle = ":", c = 'k', linewidth = linewidth)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits1 != []:
        ax[1].set_xlim(xlimits1[0], xlimits1[1])
    if ylimits != []:
        ax[1].set_ylim(ylimits[0], ylimits[1])

    ax[1].set_xticks(xticks)
    ax[1].set_yticks(yticks)

    # equation = r'$\itI_{PC} = \frac{B}{A} ln(1+AP)$'
    # ax[1].text(0.3, 0.1, equation, fontsize = fontsize+1, transform=ax[1].transAxes)
    ax[1].tick_params(axis='x', labelsize = ticksize, direction = 'in', color = 'k')
    ax[1].tick_params(axis='y', labelsize = ticksize, direction = 'in', color = 'k')
    ax[1].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
    ax[1].set_ylabel(ylabel, fontsize = fontsize, labelpad = 12, rotation = 270)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position('right')


    '''
                                                alpha
    '''

    path = "C:/QMO/Built"
    name = "2021_03_25_51_powercube_twobody_28_0.5"
    name1 = "2021_03_25_51_powercube_twobody_29_0.5"

    alpha, beta, cmap, alphaerror, betaerror, cerror, noise, yrange, xval, yval, rfi, info = loader.load_twobodyfit(path, name)
    alpha1, beta1, cmap1, alphaerror1, betaerror1, cerror1, noise1, yrange1, xval1, yval1, rfi1, info1 = loader.load_twobodyfit(path, name1)
    alpha = (alpha + alpha1) * 0.5
    noise = (noise + noise1) / 2
    yrange = (yrange + yrange1) / 2
    alpha = alpha * 0.003845

    xlabel = "$\itV_{SD}$ (V)"
    ylabel = r'$ \gamma$ / $(\beta + \alpha)$ (arb.)'
    zlabel = "Laser Power (mW)"
    fontsize = 14
    ticksize = 14
    linewidth = 5
    markersize = 15

    mapcolor_name = 'plasma' #Color of map
    colormap_name = 'viridis_r' #Color of slices

    gate = -6.5
    powerslice = [0, .01, .02, .03] # Go check cval after power arranging
    speccolor = []

    hlines = []
    # vlines = [1, 1.5]

    # fill = [0.25, 1.047, 2.4]
    # fillclr = ['deepskyblue', 'darkorange']

    smooth = 1
    noisetrim = 2
    noisetrimsdcutoff = .9


    #trims data if nonzero
    # xlimits = [.25, 2.25]
    ylimits = [-150, 3350]
    ylimits = [-.2, 18]
    ylimits = [-5.5, 12]
    xticks = [.5, 1, 1.5, 2]
    yticks = [0, 1, 5, 10]


    ############################################# DATA DO STUFF

    xlen = xval.size
    ylen = yval.size

    if noisetrim > 0:
        for i in range(xlen):
            for p in range(ylen):
                if xval[i] < noisetrimsdcutoff:
                    if noisetrim * noise[p,i] > yrange[p,i]:
                        alpha[p,i] = 0

    if smooth > 0:
        alpha = ndimage.gaussian_filter(alpha, smooth)

    powerslice = cval[:8:1]

    ############################################# DATA END STUFF

    # gateindex = np.searchsorted(yval, gate)
    y = alpha[gateindex,:]
    ax[2].plot(xval, y, linewidth = linewidth, c = 'k')

    for i in hlines:
        ax[2].axhline(i, linestyle = "dotted", c = 'k', linewidth = linewidth)
    for i in vlines:
        ax[2].axvline(i, ymin = .26, linestyle = "dotted", c = vlcolor, linewidth = linewidth, alpha = va)

    #Sets x and y range -- y range is standard set unless specified while x range is trimmed to fit data
    if xlimits != []:
        ax[2].set_xlim(xlimits[0], xlimits[1])
    if ylimits != []:
        ax[2].set_ylim(ylimits[0], ylimits[1])

    ax[2].set_xticks(xticks)
    ax[2].set_yticks(yticks)

    # for i in range(len(fill)-1):
    #     ax[2].axvspan(fill[i], fill[i+1], color = fillclr[i], alpha = filla)

    ax[2].tick_params(axis='x', labelsize = ticksize, direction = 'in')
    ax[2].tick_params(axis='y', labelsize = ticksize, direction = 'in')
    ax[2].set_xlabel(xlabel, fontsize = fontsize, labelpad = 0)
    ax[2].set_ylabel(ylabel, fontsize = fontsize, labelpad = 0)


    if save:
        pl.savefig("C:/QMO/generated plots/electron_hole/" + savename + ".png")

    if show:
        pl.show()