'''
Builds datastack from Erfu's file save system
'''
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from scipy import ndimage
from os.path import join
from os import rename
from datetime import date

from plot_tools import *
from analysis_tools import *

path = "E:/Data/Raw/Lui_PL_Data/MW32_SD"
name = "M14_"
newname = "MW32_150K_plstack"
savepath = "C:/Jed/Built"

vsd = np.linspace(-5,5,101)
vg = np.linspace(-5,5,101)
d = np.loadtxt(join(path, name + str(0).zfill(3) + ".csv"), delimiter=",")
wval = d[:,0]

#the number of times wval is wrong (hopefully none but this is to check)
ecount = 0

data = np.zeros((vsd.size, vg.size, wval.size))

#Loop through every file and load it into memory
for p in range(vg.size):
    for i in range(vsd.size):
        count = (p*vsd.size) + i
        fname = name + str(count).zfill(3) + ".csv"
        d = np.loadtxt(join(path, fname), delimiter=",")
        #d[:,0] is wval and d[:,1] is the data

        if d[12,0] != wval[12]:
            ecount += 1
        if i == 0:
            print(fname + "   " + str(ecount))
        data[i,p,:] = d[:,1]


rf = np.zeros(d.shape)
power = np.zeros(d.shape)
info = np.zeros((1))
savename = join(savepath, newname)
np.savez(savename, d = data, xval = vsd, yval = vg, cval = wval, rfi=rf, pow = power, inf = info)