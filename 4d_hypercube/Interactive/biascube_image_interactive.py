'''
A script to actively explore a hypercube
Displays a map of the source and gate on one figure as well as correctable colorbars
Second figure shows pc image and pc conductance image
'''
from scipy import ndimage
from scipy.ndimage.filters import gaussian_filter

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as pl

from datetime import date
import sys
sys.path.append("C:/Users/Jedki/Documents/analCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
# from builder import build_map

path = "C:/Users/jedki/QMOdata/Built"
# name = "2020_08_12_biascube"
name = "2020_09_29_biascube_noa"
name = "2020_09_03_biascube"
name = "2021_10_29_powermap_noa"     
xval, yval, cval, zval, data, rf = loader.load_hypercube(path, name)

pcicolor = "plasma"
conductioncolor = "magma"
biasmapcolor = 'plasma'

global recolor
global cmin
global cmax
global zmin
global zmax
recolor = False
contours = True
zmin = 0
zmax = 0
dmin = 0
dmax = 0

# Build axes
axs, axcbs, figs = makeaxes(2, cbar = True)
fi = pl.figure(figsize=figs)
axs1,axcbs1, figs = makeaxes(1, cbar = True)
fig = pl.figure(figsize = figs)

ax1 = fi.add_axes(axs[0])
ax2 = fi.add_axes(axs[1])
axcb = fi.add_axes(axcbs[0])

a1 = fig.add_axes(axs1[0])
axcb1 = fig.add_axes(axcbs1[0])

#Build colormaps and set initial values

Vsd = cval
Vg = zval
zmin, zmax, cmin, cmax = 0,0,0,0

global gateindex
gateindex = int(Vg.size / 2)
global sdindex
sdindex = int(Vsd.size / 2)

derstep = cval[0] - cval[1]

cdata = np.gradient(data, derstep, axis = 2)

#Initial map states
meanmargin = 15
cmargin = 10
glen = Vg.size
sdlen = Vsd.size
xlen = xval.size
ylen = yval.size
sdstep = np.abs((Vsd[0]-Vsd[-1])/sdlen)

zd = data[:,:,sdindex, gateindex]
cd = cdata[:,:,sdindex, gateindex]
biasmap = np.zeros(data[0,0,:,:].shape)
if contours:
    rfd = rf[:,:,sdindex, gateindex]

for i in range(glen):
    for n in range(sdlen):
        biasmap[n,i] = np.mean(data[meanmargin:ylen-meanmargin,meanmargin:xlen-meanmargin,n,i])
biasmap = np.transpose(biasmap)

#colormap colors
if recolor:
    pcmap, pcNorm = color_memap( pcicolor,zd[cmargin:ylen-cmargin, cmargin:xlen-cmargin],dmin = zmin, dmax = zmax)
    ccmap, ccNorm = color_memap( conductioncolor, cd[cmargin:ylen-cmargin, cmargin:xlen-cmargin], dmin = cmin, dmax = cmax)
else:
    pcmap, pcNorm = color_memap( pcicolor,zd, dmin = zmin, dmax = zmax)
    ccmap, ccNorm = color_memap( conductioncolor, cd, dmin = cmin, dmax = cmax)
bcmap, bcNorm = color_memap( biasmapcolor, biasmap, dmin = 0, dmax = 0)

xt = get_extent(xval, yval)
bt = get_extent(Vsd, Vg)

cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=pcmap, norm=pcNorm, orientation='vertical')
axcb.set_title(r'$I (nA)$')

ax1.imshow(zd, cmap=pcmap, norm=pcNorm, extent = xt, interpolation='none', origin = 'lower')
ax2.imshow(cd, cmap=ccmap, norm=ccNorm, extent = xt, interpolation='none', origin = 'lower')
if contours:
    ax1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
    ax2.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
# ax1.set_title(r"$V_{g}$ = " + str(round(zVg,1)) + r"V, $V_{sd}$ = " + str(round(zVsd, 1)) + 'mV')
# point1 = ax1.plot([ixx], [ixy], 'ro')

a1.imshow(biasmap, extent=bt, cmap=bcmap, norm=bcNorm, interpolation='none', aspect='auto', origin = 'lower')
a1.scatter(Vsd[sdindex],Vg[gateindex], color = 'k', marker = 's')
fig.suptitle("Source = " + str(Vsd[sdindex]) + "V  Gate = " + str(Vg[gateindex]))

ax2.set_xlabel(r"V$_{sd}$ (mV)")
ax2.set_ylabel(r"V$_{g}$ (V)")

def onClick(event):
    global sdindex
    global gateindex
    global recolor
    global cmin
    global cmax
    global zmin
    global zmax
    if event.inaxes == a1:
        ix = int(event.xdata)
        iy = int(event.ydata)
        print(event.xdata)
        print(event.ydata)        
        sdindex = np.searchsorted(Vsd, ix)
        gateindex = np.searchsorted(Vg, iy)
        if sdindex == sdlen:
            sdindex = sdindex - 1
        if gateindex == glen:
            gateindex = gateindex - 1
        zd = data[:,:,sdindex, gateindex]
        cd = cdata[:,:,sdindex, gateindex]
        if recolor:
            pcmap, pcNorm = color_memap( pcicolor,zd[cmargin:ylen-cmargin, cmargin:xlen-cmargin])
            ccmap, ccNorm = color_memap( conductioncolor, cd[cmargin:ylen-cmargin, cmargin:xlen-cmargin])
        else:
            pcmap, pcNorm = color_memap( pcicolor,zd, dmin = zmin, dmax = zmax)
            ccmap, ccNorm = color_memap( conductioncolor, cd, dmin = cmin, dmax = cmax)
        ax1.cla()
        ax2.cla()
        a1.cla()
        axcb.cla()
        ax1.imshow(zd, cmap=pcmap, norm=pcNorm, extent = xt, origin = 'lower')
        ax2.imshow(cd, cmap=ccmap, norm=ccNorm, extent = xt, origin = 'lower')
        a1.imshow(biasmap, extent=bt, cmap=bcmap, norm=bcNorm, aspect='auto', origin = 'lower')
        a1.scatter(Vsd[sdindex],Vg[gateindex], color = 'k', marker = 's')
        fig.suptitle("Source = " + str(Vsd[sdindex]) + "V  Gate = " + str(np.around(Vg[gateindex], 3)))
        cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=pcmap, norm=pcNorm, orientation='vertical')
        fi.canvas.draw()
        fig.canvas.draw()
        fig.canvas.flush_events()

def onKey(event):
    global sdindex
    global gateindex
    global recolor
    global cmin
    global cmax
    global zmin
    global zmax
    if event.key == "left":
        if sdindex > 0:
            sdindex = sdindex - 1
    if event.key == "right":
        if sdindex < sdlen-1:
            sdindex = sdindex + 1
    if event.key == "up":
        if gateindex < glen - 1:
            gateindex = gateindex + 1
    if event.key == "down":
        if gateindex > 0:
            gateindex = gateindex - 1

    zd = data[:,:,sdindex, gateindex]
    cd = cdata[:,:,sdindex, gateindex]
    rfd = rf[:,:,sdindex, gateindex]
    if recolor:
        pcmap, pcNorm = color_memap( pcicolor,zd[cmargin:ylen-cmargin, cmargin:xlen-cmargin])
        ccmap, ccNorm = color_memap( conductioncolor, cd[cmargin:ylen-cmargin, cmargin:xlen-cmargin])
    else:
        pcmap = pl.get_cmap(pcicolor)
        pcNorm = mpl.colors.Normalize(vmin = zmin, vmax = zmax)
        ccmap = pl.get_cmap(conductioncolor)
        ccNorm = mpl.colors.Normalize(vmin = cmin, vmax = cmax)
    ax1.cla()
    ax2.cla()
    a1.cla()
    axcb.cla()
    ax1.imshow(zd, cmap=pcmap, norm=pcNorm, extent = xt, origin = 'lower')
    ax2.imshow(cd, cmap=ccmap, norm=ccNorm, extent = xt, origin = 'lower')
    if contours:
        ax1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5,origin = 'lower' )#, aspect = 'auto')# origin = 'lower')
        ax2.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5,origin = 'lower' )#, aspect = 'auto')# )#)
    a1.imshow(biasmap, extent=bt, cmap=bcmap, norm=bcNorm, aspect='auto', origin = 'lower')
    a1.scatter(Vsd[sdindex],Vg[gateindex], color = 'k', marker = 's')
    fig.suptitle("Source = " + str(Vsd[sdindex]) + "V  Gate = " + str(np.around(Vg[gateindex], 3)))
    cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=pcmap, norm=pcNorm, orientation='vertical')
    fi.canvas.draw()
    fig.canvas.draw()
    fig.canvas.flush_events()



# end onClick

fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('key_press_event', onKey)


fig.show()
pl.show()
