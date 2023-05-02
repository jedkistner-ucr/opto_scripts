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

path = "C:/QMO/Raw Data/Do3_deltaV"
newname = "ao3_deltav"
savepath = "C:/QMO/Built"

#No prior knowledge
yval = []
wval = []
data = []

#pull every file in the path
f = listdir(path)
ecount = 0
tcount = 0

for file in f:

    if data == []:
        d = np.loadtxt(join(path, file), skiprows=1)
        wval = d[:,0]

    name = file.split("_")
    v1 = name[1].split("=")[1]
    v1 = float(v1[:len(v1)-5])


    yval.append(v1)

    d = np.loadtxt(join(path, file), skiprows=1)
    data.append(d[:,1])

    if wval[12] != d[12,0]:
        ecount += 1
        print(ecount)

    if tcount % 10 == 0:
        print(tcount)
    tcount += 1

#data and yval and now filled, we just need to sort them
yval_sort = sorted(yval)
xy = sorted(zip(yval, data))
data = [x for y, x in xy]

yval = np.asarray(yval)
print(yval_sort)
data = np.asarray(data)

savename = join(savepath, newname)
np.savez(savename, d = data, xval = wval, yval = yval_sort)