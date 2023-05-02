'''
A script to actively explore a cube
Choose a point on the map and view the z axis liencut
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
name = "2021_10_21_14_powercube" 
name = "2022_01_26_71_powercube"
name = "2022_02_02_48_powercube"
name = "2022_02_02_51_powercube"
name = "2022_02_03_32_powercube"
name = "2022_02_03_12_powercube"
xval, yval, cval, data, rf, power, info = loader.load_cube(path, name)

pcicolor = "plasma"
conductioncolor = "magma"
biasmapcolor = 'plasma'

gate = 0
sd = 0

usepower = True
zerooffset = True

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
contours = False
smooth = False
smoothval = 1
smoothconduction = False
smoothcval = 1
conductionmap = False
invertconduction = True
zmin = 0
zmax = .02
dmin = 0
dmax = .3

clen = cval[:].size
xlen = xval[:].size
ylen = yval[:].size

# Build axes
axs, figs = makeaxes(1, cbar = False)
fi = pl.figure(figsize=figs)
axs1,axcbs1, figs = makeaxes(1, cbar = True)
fig = pl.figure(figsize = figs)

ax1 = fi.add_axes(axs[0])
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

if usepower:
    for i in range(cval.size):
        cval[i] = np.mean(power[:,:,i])
if zerooffset:
    offset = np.mean(data[:,:,0])
    data = data - offset


#bias index for clicking map
cindex = clen - 1# np.searchsorted(cval, sd)

#image and contours for clickable map
pcimage = data[:,:,cindex]
if contours:
    rfd = rf[:,:,cindex]

#starting bias map on display figure
d = data[yindex,xindex,:]

dmap, dNorm = color_memap( pcicolor , pcimage )
pcmap, pcNorm = color_memap( pcicolor , pcimage )

xt = get_extent(xval, yval, inverty= True)

#Set initial plot values
#three plots on clickable figure
a1.imshow(pcimage, extent=xt, cmap=pcmap, norm=dNorm, interpolation='none', aspect='auto', origin = 'lower')
if contours:
    a1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5, aspect = 'auto', origin = 'lower')
a1.scatter(xval[xindex],yval[yindex], color = 'k', marker = 's')

x = cval
y = d

ax1.plot(x, y, linewidth = 1, c = 'k')
# ax3.plot(x, yc, linewidth = 1, c = clr[i])

ax1.set_xlim(cval[0], cval[-1])
# ax3.set_xlim(cval[0], cval[-1])

fig.suptitle("Photocurrent Map")

ax1.set_xlabel(r"V$_{sd}$ (mV)")
ax1.set_ylabel(r"V$_{g}$ (V)")

def onKey(event):
    global xindex
    global yindex
    global gateindex
    global cindex
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

    zd = data[yindex,xindex,:]

    rfd = rf[:,:,cindex]

    ax1.cla()
    a1.cla()
    # ax3.cla()

    x = cval
    y = zd

    ax1.plot(x, y, linewidth = 1, c = 'k')
    # ax3.plot(x, yc, linewidth = 1, c = clr[i])

    a1.imshow(pcimage, extent=xt, cmap=pcmap, norm=pcNorm, interpolation='none', aspect='auto', origin = 'lower')
    # if contours:
    #     a1.contour( rfd, cmap = pl.get_cmap('Greys'), linewidths = 1,  extent = xt, levels = 5,origin = 'lower' )#, aspect = 'auto')# origin = 'lower')
    a1.scatter(xval[xindex],yval[yindex], color = 'k', marker = 's')
    
    ax1.set_xlim(cval[0], cval[-1])
    # ax3.set_xlim(cval[0], cval[-1])

    # fig.suptitle("Source = " + str(Vsd[cindex]) + "V  Gate = " + str(np.around(Vg[gateindex], 3)))
    fi.canvas.draw()
    fig.canvas.draw()
    fig.canvas.flush_events()



# end onClick

# fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('key_press_event', onKey)


fig.show()
pl.show()
