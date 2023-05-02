'''
Builds datastack from Ao's file save system (really erfu's but this one is saved differently)
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from scipy import ndimage
from os.path import join
from os import rename, listdir
from datetime import date
import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *

path = "C:/QMO/Raw Data/M9"
newname = "erfu01_tg_bg"
savepath = "C:/QMO/Built"

#Prior knowledge
xval = []
yval = []
wval = []
data = []

#pull every file in the path
f = listdir(path)
ecount = 0
tcount = 0

#generate xval and yval from filenames
for file in f:
    name = file.split("_")
    v1 = name[1].split("=")[1]
    v2 = name[2].split("=")[1]
    v1 = float(v1[:len(v1)-1])
    v2 = float(v2[:len(v2)-5])
    # print(str(v1) + "  " + str(v2))

    xval.append(v1)
    yval.append(v2)

xval = np.sort(np.unique(np.asarray(xval)))
yval = np.sort(np.unique(np.asarray(yval)))

#loop through files again and load data
for file in f:

    if data == []:
        d = np.loadtxt(join(path, file), skiprows=1, delimiter = '\t')
        wval = d[:,0]
        data = np.zeros((xval.size, yval.size, wval.size))

    name = file.split("_")
    v1 = name[1].split("=")[1]
    v2 = name[2].split("=")[1]
    v1 = float(v1[:len(v1)-1])
    v2 = float(v2[:len(v2)-5])

    #Find the appropriate index for the file
    v1index = np.searchsorted(xval, v1)
    v2index = np.searchsorted(yval, v2)

    d = np.loadtxt(join(path, file), skiprows=1, delimiter = '\t')
    data[v1index, v2index,:] = d[:,1]

    if wval[12] != d[12,0]:
        ecount += 1
        print(ecount)

    if tcount % 10 == 0:
        print(tcount)
    tcount += 1

savename = join(savepath, newname)
np.savez(savename, d = data, xval = xval, yval = yval, cval = wval)