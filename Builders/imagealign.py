'''
takes a set of rfi/pci images and alignes them
uses a brute force technique where it finds the minimum difference between each image and the starting image
outputs:
data[x, y, cube , repeat]  rf[x, y, cube , repeat]
xval[] yval[] single dimension array taken from first scan as everything else is aligned to that
cval[] zval[] single dimension cube axis and repeat axis
i decided not to fuck with power and info
'''

import numpy as np
from os.path import join
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
from scipy import ndimage
from plot_tools import *
from analysis_tools import *

############ Run parameters

path = "C:/QMO/Raw Data/2023_01/2023_01_31" #data location where scans are stored
newpath = "C:/QMO/Built"    #new data file save location
newname = "2023_01_30_biasmap_noa"                #new data file save name

#The name of the repeating variable in the log file -- this will become zval
repeatname = "Backgate Start"
# repeatname = "Source/Drain Start"

# Build a list of the data files that you want to allign
#03_25_68 to 03_26_17
# names = ['CPS_2021_03_25_68']
# for i in range(8):
#     num = 69 + i
#     names.append('CPS_2021_03_25_' + str(num))
# for i in range(18):
#     num = 0 + i
#     names.append('CPS_2021_03_26_' + str(num))

# Build a list of the data files that you want to allign
#03_24_74 to 03_25_23
# names = ['CPS_2021_03_24_74']
# for i in range(4):
#     num = 75 + i
#     names.append('CPS_2021_03_24_' + str(num))
# for i in range(24):
#     num = 0 + i
#     names.append('CPS_2021_03_25_' + str(num))

#65-8
names = []
for i in range(9):
    num = 0 + i
    names.append('CPS_2023_01_31_' + str(num))


#Build a list of the data files that you want to allign
# names = ['CPS_2021_10_28_44']
# for i in range(6):
#     num = 45 + i
#     names.append('CPS_2021_10_28_' + str(num))
# for i in range(14):
#     num = 0 + i
#     names.append('CPS_2021_10_29_' + str(num))

# Build a list of the data files that you want to allign
# names = ['CPS_2022_01_27_7']
# for i in range(8):
#     num = 8 + i
#     names.append('CPS_2022_01_27_' + str(num))
# for i in range(17):
#     num = 0 + i
#     names.append('CPS_2022_01_28_' + str(num))
# for i in range(7):
#     num = 0 + i
#     names.append('CPS_2022_01_29_' + str(num))

    # Build a list of the data files that you want to allign
# names = ['CPS_2022_01_26_33']
# for i in range(10):
#     num = 34 + i
#     names.append('CPS_2022_01_26_' + str(num))



#Alignment parameters
#The margin in pixels to do the alignment over: margin = 10 will check an 11 x 11 grid
xmargin = 0
ymargin = 0
userfi = True #use rf or pci for alignment
smooth = True #Smooths the image
smoothval = 0
cutlowvalues = False #Values below a threshold are set to zero -- only works for rf right now
cutoff = -.003
tiltcorrect = False #Tilts the subtraction matrix to correct for abnormalities
zerobackground = False
backgroundspot = [42,0]
bgmargin = -.5

#index of image to align other images to (the primary image) -1 will be set to the end of the stack
alignmentindex = [-1, 0]

#option to average power data files for cval instead of recorded logfile values
usepower = False
usepvalhack = False
#Pixel values to do the alignment over, effectively crops the images to avoide edge effects
#if both min and max are zero it will use the full image
xmin = 0
xmax = 0
ymin = 0
ymax = 0

############ End run parameters

#compiles images into a single array
#First pull logfile to get xval and yval
logfile = path + "/" + names[0] + "_log.log"
f = open(logfile, "r")
info = {}
info["scan"] = names[0]

for line in f:
        s = line.split(":")
        info[s[0]] = s[1]
f.close()

xval = np.linspace( float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace( float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"]))
cval = np.linspace( float(info["Cube Axis Start"]), float(info["Cube Axis End"]), int(info["Cube Scan Number"]))
zval = []

if usepvalhack:
    cval = getpvalhack(cval)

#Zval is pulled from successive log files and rf and data are pulled from successive .npy files
datalist = []
rflist = []
powlist = []
for i in range(len(names)):
    #Pull pci and rfi files
    d = np.load(join(path, names[i] + "_" + "pci" + ".npy"))
    r = np.load(join(path, names[i] + "_" + "rfi" + ".npy"))


    p = np.load(join(path, names[i] + "_" + "pow" + ".npy"))
    
    datalist.append(d)
    rflist.append(r)
    powlist.append(p)

    logfile = path + "/" + names[i] + "_log.log"
    f = open(logfile, "r")
    for line in f:
        s = line.split(":")
        if s[0] == repeatname:
            zval.append(float(s[1]))


if usepower:
    p = np.load(join(path, names[0] + "_" + "pow" + ".npy"))
    for i in range(cval.size):
        cval[i] = np.mean(p[:,:,i])


#Data is stacked
data = np.stack(datalist, axis = -1)
rf = np.stack(rflist, axis = -1)
power = np.stack(powlist, axis = -1)
zval = np.asanyarray(zval)

data = data * float(info["Pre-Amp Gain"])

#Sets alignment indices if any are -1
if alignmentindex[0] == -1:
    alignmentindex[0] = cval.size - 1
if alignmentindex[1] == -1:
    alignmentindex[1] = zval.size - 1

#Empty arrays to fill with the corrected data
aldata = np.zeros(data.shape)
alrf = np.zeros(data.shape)
alpower = np.zeros(data.shape)

offsets = []    #unused, for debug or curiosity, a list of the offsets


#Sets values below the cutoff to zero
if cutlowvalues:
    for n in range(rf[0,0,:,0].size):
        for m in range(rf[0,0,0,:].size):
            for x in range(rf[0,:,0,0].size):
                for y in range(rf[:,0,0,0].size):
                    if rf[y, x, n, m] < cutoff:
                        rf[y, x, n, m] = 0

#Trims the pixel area that the alignment will be done over
if xmin == 0 and xmax == 0:
    x0 = 0
    x1 = data[0,:,0,0].size
else:
    x0 = xmin
    x1 = xmax
if ymin == 0 and ymax == 0:
    y0 = 0
    y1 = data[:,0,0,0].size
else:
    y0 = ymin
    y1 = ymax
errorlist = []

#Start alignment, looping though every image in the set
primary = np.zeros(data[y0:y1,x0:x1,0,0].shape)
secondary = np.zeros(data[y0:y1,x0:x1,0,0].shape)
trimprim = np.zeros(data[y0:y1,x0:x1,0,0].shape)
nlen = data[y0:y1,0,0,0].size   #lengths to save allocation time later
mlen = data[0,x0:x1,0,0].size

#Sets primary image
if userfi:
    if smooth:
        primary = applyGauss(rf[y0:y1,x0:x1,alignmentindex[0],alignmentindex[1]], smoothval)
    else:
        primary = rf[y0:y1,x0:x1,alignmentindex[0],alignmentindex[1]]
else:
    primary[:,:] = data[y0:y1,x0:x1,alignmentindex[0],alignmentindex[1]]
    if zerobackground:  #zeros the background
        bgvalue = 0
        for bi in range(2):
            for bn in range(2):
                bgvalue = bgvalue + primary[backgroundspot[1] + bi, backgroundspot[0] + bn]
                # print()
        bgvalue = bgvalue / 4
        
        for x in range(primary[0,:].size):
            for y in range(primary[:,0].size):
                if primary[y, x] > bgvalue + bgmargin:
                    primary[y,x] = 0
    
# pl.figure()
# pl.imshow(primary)
# pl.scatter(backgroundspot[0], backgroundspot[1])
# print(primary[backgroundspot[1], backgroundspot[0]])
# print(bgvalue)
# pl.show()

for i in range(data.shape[2]):
    for p in range(data.shape[3]):

        if i == alignmentindex[0] and p == alignmentindex[1]:
            alrf[:,:,i,p] = rf[:,:,i,p]
            aldata[:,:,i,p] = data[:,:,i,p]
            alpower[:,:,i,p] = power[:,:,i,p]
        else:   #Otherwise builds matrix of image subtraction values
            diff = np.ones((2*ymargin+1 , 2*xmargin+1))
            diff = diff * 100
            testn = -1
            

            if userfi:
                if smooth:
                    secondary = applyGauss( rf[y0:y1,x0:x1,i,p], smoothval )
                else:
                    secondary = rf[y0:y1,x0:x1,i,p]
            else:
                secondary[:,:] = data[y0:y1,x0:x1,i,p]

                if zerobackground:
                    bgvalue = 0
                    for bi in range(2):
                        for bn in range(2):
                            bgvalue = bgvalue + secondary[backgroundspot[1] + bi, backgroundspot[0] + bn]
                            # print()
                    bgvalue = bgvalue / 4
        
                    for x in range(secondary[0,:].size):
                        for y in range(secondary[:,0].size):
                            if secondary[y, x] > bgvalue + bgmargin:
                                secondary[y,x] = 0

            # pl.figure()
            # pl.imshow(secondary)
                    # pl.show()

            #For every position in the offset matrix
            for n in range(-ymargin, ymargin +1,  1):
                for m in range(-xmargin, xmargin +1,  1):


                    trimprim[:,:] = primary[:,:]
                    shiftedrf = ndimage.shift(secondary, [n, m])

                    #Trims the primary image so that not considered areas are zero
                    if n > 0:
                        for nn in range(n ):
                            trimprim[nn ,:] = np.zeros(trimprim[0,:].shape)
                    elif n < 0:
                        for nn in range(n,0):
                            trimprim[nlen + nn  ,:] = np.zeros(trimprim[0,:].shape)
                            # print()

                    if m > 0:
                        for mm in range(m):
                            trimprim[:, mm] = np.zeros(trimprim[:,0].shape)
                    elif m < 0:
                        for mm in range(m, 0):
                            trimprim[: , mlen + mm] = np.zeros(trimprim[:,0].shape)

                    #Builds array of error
                    subtracted = np.abs(trimprim - shiftedrf)
                    sumsub = np.sum(subtracted)
                    totpixs =  (nlen - np.abs(n))  *  (mlen - np.abs(m))   
                    meansub = sumsub / totpixs
                    diff[n + ymargin, m + xmargin] = meansub

                    # if n > 2:
                    #     pl.figure()
                    #     pl.imshow(primary)
                    #     # pl.show()
                    #     pl.figure()
                    #     pl.imshow(secondary)
                    #     pl.figure()
                    #     pl.imshow(trimprim)
                    #     # pl.show()
                    #     pl.figure()
                    #     pl.imshow(shiftedrf)
                    #     pl.show()


            #Gets minimum index from error array
            minindex = np.unravel_index(diff.argmin(), diff.shape)
            shiftval = [(minindex[0] - ymargin), (minindex[1] - xmargin)]
            offsets.append(shiftval)
            # errorlist.append(diff)
            # if p > -1:
            #     print(" ")
            #     print(shiftval)
            # pl.figure()
            # pl.imshow(diff)
            # pl.scatter(minindex[1], minindex[0])
            # pl.show()

            #Fills alrf and aldata with shifted data
            alrf[:,:,i,p] = ndimage.shift(rf[:,:,i,p], shiftval)
            aldata[:,:,i,p] = ndimage.shift(data[:,:,i,p], shiftval)
            alpower[:,:,i,p] = ndimage.shift(power[:,:,i,p], shiftval)


            # if np.abs(shiftval[0]) > 4 or np.abs(shiftval[1]) > 4: 
            # if p > 0:
            #     print("")
            #     print(minindex)
            #     print(shiftval)

            #     pl.figure()
            #     pl.imshow(ndimage.shift(data[:,:,i,p], shiftval))
            #     pl.contour(primary)
            #     pl.show()


            # print(p)
    print(i)    #running tracker so you know its working, takes about 20 seconds per repeat cube

#saves
savename = join(newpath, newname)
np.savez(savename, data = aldata, xval = xval, yval = yval, cval = cval, zval = zval, rfi = alrf, pow = alpower)