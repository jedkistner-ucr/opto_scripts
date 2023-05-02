'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

from matplotlib import use
import numpy as np
from os.path import join
import sys
from analysis_tools import *


usepower = False
year = 2021
month = 11
day = 3
run = 2

path = "E:/Data/Raw/%i/%i_%.2i/%i_%.2i_%.2i/" % (year, year, month, year, month, day)
newpath = "C:/Jed/Built"    #new data file save location
name = "cps_%i_%.2i_%.2i_%i" % (year, month, day, run)
date = "%i_%.2i_%.2i_%i" % (year, month, day, run)
newname = date + "_sdcube" 

print(newname)

datapath = join(path, name + "_pci.npy")
logpath = join(path, name + "_log.log")
powpath = join(path, name + "_pow.npy")
rfpath = join(path, name + "_rfi.npy")

data = np.load(datapath)
powdata = np.load(powpath)
rfdata = np.load(rfpath)

info = {}
info["Scan"] = date

f = open(logpath, "r")

for line in f:
        s = line.split(":")
        info[s[0]] = s[1]
        if s[0] == "DATA_dir_path":
                break
f.close()

xval = np.linspace( float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace( float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"])) 

if usepower:
    pval = []
    for i in range(powdata[0,0,:].size):
        pval.append(np.mean(powdata[:,:,i]))

    pval = np.asanyarray(pval)
else:
    pval = np.linspace( float(info["Cube Axis Start"]), float(info["Cube Axis End"]), int(info["Cube Scan Number"]))

data = data * float(info["Pre-Amp Gain"])

xlen = xval.size
ylen = yval.size
plen = pval.size

inf = np.empty((2,len(info)), dtype='object')
i = 0
for k in info:
    inf[0][i] = k
    inf[1][i] = info[k]
    i = i+1

savename = join(newpath, newname)
np.savez(savename, d = data, xval = xval, yval = yval, cval = pval, rfi = rfdata, pow = powdata, inf = inf)