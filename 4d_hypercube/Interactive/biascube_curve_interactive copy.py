'''
A script to actively explore a hypercube
Shows an image of the photocurrent on one figure
Second figure has conduction map, transport, and conduction curve
Space "saves" a curve and keeps it in the bacground as gray curves
Backspace clears saved curves
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
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
# from builder import build_map

path = "C:/QMO/Built"
# name = "2020_08_12_biascube"
name = "2021_03_17_powercube_noa"
name = "2022_01_29_powermap_noa"
xval, yval, cval, zval, data, rf, power1 = loader.load_hypercube(path, name)
# zval = np.flip(zval, 0)
# data = np.flip(data, 3)
# cval = cval[1:-1]
# data = data[:,:,1:-1,:]

pcicolor = "seismic"
conductioncolor = "magma"
biasmapcolor = 'plasma'

gate = 0
sd = 0

global rescale
global recolor
global cmin
global cmax
global zmin
global zmax
global ymin
global ymax
recolor = False
rescale = False
contours = True
smooth = True
smoothval = 1
smoothconduction = True
smoothcval = 1
conductionmap = False
invertconduction = True
zmin = -.08
zmax = .08
dmin = -.06
dmax = .06

clen = cval[:].size
zlen = zval[:].size
xlen = xval[:].size
ylen = yval[:].size

# Build axes
axs, axcbs, figs = makeaxes(3, cbar = True)
fi = pl.figure(figsize=figs)
axs1,axcbs1, figs = makeaxes(1, cbar = True)
fig = pl.figure(figsize = figs)

ax1 = fi.add_axes(axs[0])
ax2 = fi.add_axes(axs[1])
ax3 = fi.add_axes(axs[2])
axcb = fi.add_axes(axcbs[0])

a1 = fig.add_axes(axs1[0])
axcb1 = fig.add_axes(axcbs1[0])

#Data manipulation

#Build colormaps and set initial values
zmin, zmax, cmin, cmax, ymin, ymax = 0,0,0,0,0,0

#initial spot for x and y
global yindex
yindex = int(ylen / 2)
global xindex
xindex = int(xlen / 2)

#generates conduction map
derstep = np.abs(zval[0] - zval[1])

if smooth:
    for i in range(xlen):
        for p in range(ylen):
            data[p,i,:,:] = applyGauss( data[p,i,:,:], smoothval)

cdata = np.gradient( data, derstep, axis = 3)
if invertconduction:
    cdata = cdata * -1
#Initial map states ...

#bias index for clicking map
sdindex = 25 #np.searchsorted(cval, sd)
gateindex = 15 #np.searchsorted(zval, gate)

#image and contours for clickable map
pcimage = data[:,:,sdindex, gateindex]
if contours:
    rfd = rf[:,:,sdindex, gateindex]

#starting bias map on display figure
d = data[yindex,xindex,:, :]
cd = cdata[yindex,xindex,:, :]

if smoothconduction:
    cd = applyGauss(cd, smoothcval)

#get colors
if recolor:
    cmargin = 2
    dmap, dNorm = color_memap( pcicolor,d[cmargin:clen-cmargin, cmargin:zlen-cmargin] )
    cdmap, cdNorm = color_memap( conductioncolor, cd[cmargin:clen-cmargin, cmargin:zlen-cmargin] )
else:
    dmap, dNorm = color_memap( pcicolor,d, dmin = zmin, dmax = zmax)
    cdmap, cdNorm = color_memap( conductioncolor, cd, dmin = cmin, dmax = cmax)

pcmap, pcNorm = color_memap( pcicolor,pcimage )

xt = get_extent(xval, yval)
bt = get_extent(cval, zval)

if conductionmap:
    cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=cdmap, norm=cdNorm, orientation='vertical')
else:
    cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=dmap, norm=dNorm, orientation='vertical')

clr = color_meline('plasma', zlen)

axcb.set_title(r'$I (nA)$')

#Set initial plot values
#three plots on clickable figure
a1.imshow(pcimage, extent=xt, cmap=pcmap, norm=dNorm, interpolation='none', aspect='auto', origin = 'lower')
if contours:
    a1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
a1.scatter(xval[xindex],yval[yindex], color = 'k', marker = 's')

#plots on display figure
if conductionmap:
    ax1.imshow(cd, cmap=cdmap, norm=cdNorm, extent = bt,aspect = 'auto', interpolation='none', origin = 'lower')
else:
    ax1.imshow(d, cmap=dmap, norm=dNorm, extent = bt, aspect = 'auto', interpolation='none', origin = 'lower')

for i in range(zlen):
    x = cval
    y = d[ :, i]
    yc = cd[ :, i]

    ax2.plot(x, y, linewidth = 1, c = clr[i])
    ax3.plot(x, yc, linewidth = 1, c = clr[i])

ax2.set_xlim(cval[0], cval[-1])
ax3.set_xlim(cval[0], cval[-1])

fig.suptitle("Photocurrent Image")

ax2.set_xlabel(r"V$_{sd}$ (mV)")
ax2.set_ylabel(r"V$_{g}$ (V)")



def onClick(event):
    global xindex
    global yindex
    global recolor
    global rescale
    global cmin
    global cmax
    global zmin
    global zmax
    if event.inaxes == a1:
        ix = int(event.xdata)
        iy = int(event.ydata)
        print(event.xdata)
        print(event.ydata)        
        # sdindex = np.searchsorted(Vsd, ix)
        # gateindex = np.searchsorted(Vg, iy)
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
    global xindex
    global yindex
    global gateindex
    global sdindex
    global recolor
    global rescale
    global cmin
    global cmax
    global zmin
    global zmax
    if event.key == "left":
        if xindex > 0:
            xindex = xindex - 1
    if event.key == "right":
        if xindex < xlen -  1:
            xindex = xindex + 1
    if event.key == "down":
        if yindex < ylen - 1:
            yindex = yindex + 1
    if event.key == "up":
        if yindex > 0:
            yindex = yindex - 1

    zd = cdata[yindex,xindex,:, :]
    cd = cdata[yindex,xindex,:, :]
    print(sdindex)
    if smoothconduction:
        cd = applyGauss(cd, smoothcval)

    rfd = rf[:,:,sdindex, gateindex]

    if recolor:
        if conductionmap:
            dmap, dNorm = color_memap( pcicolor,cd )
        else:
            dmap, dNorm = color_memap( pcicolor, zd )
    else:
        dmap = pl.get_cmap(pcicolor)
        dNorm = mpl.colors.Normalize(vmin = zmin, vmax = zmax)
        dmap = pl.get_cmap(conductioncolor)
        dNorm = mpl.colors.Normalize(vmin = cmin, vmax = cmax)

    ax1.cla()
    ax2.cla()
    ax3.cla()
    a1.cla()
    axcb.cla()

    if conductionmap:
        ax1.imshow(zd, cmap=dmap, norm=dNorm, extent = bt, aspect = 'auto', origin = 'lower')
    else:
        ax1.imshow(cd, cmap=dmap, norm=dNorm, extent = bt, aspect = 'auto', origin = 'lower')

    for i in range(zlen):
        x = cval
        y = zd[ :, i]
        ax2.plot(x, y, linewidth = 1, c = clr[i])

    for i in range(clen):
        x1 = zval
        yc = cd[ i, :]
        ax3.plot(x1, yc, linewidth = 1, c = clr[i])

    a1.imshow(pcimage, extent=xt, cmap=pcmap, norm=pcNorm, interpolation='none', aspect='auto', origin = 'lower')
    if contours:
        a1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5,origin = 'lower' )#, aspect = 'auto')# origin = 'lower')
    a1.scatter(xval[xindex],yval[yindex], color = 'k', marker = 's')
    
    ax2.set_xlim(cval[0], cval[-1])
    ax3.set_xlim(zval[0], zval[-1])

    # fig.suptitle("Source = " + str(Vsd[sdindex]) + "V  Gate = " + str(np.around(Vg[gateindex], 3)))
    cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=pcmap, norm=pcNorm, orientation='vertical')
    fi.canvas.draw()
    fig.canvas.draw()
    fig.canvas.flush_events()



# end onClick

# fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('key_press_event', onKey)


fig.show()
pl.show()
