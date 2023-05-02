'''
Takes a just finished run, builds dataset, creates a map and linecuts then combines it with a second dataset
'''

import numpy as np
from os.path import join
import sys
sys.path.append("C:/Users/Jedki/Documents/analCode/Toolbox")
from analysis_tools import *

path = "C:/Users/Jedki/QMOdata/Raw Data/"
newpath = "C:/Users/Jedki/QMOdata/Built"    #new data file save location


dates = ["2021_10_28_20", "2021_10_28_21", "2021_10_28_22"]
prefix = "cps_"

data = []
powdata = []
rfdata = []
cval = []


for date in dates:

    newname = dates[0] + "_powercube" 
    usepower = False

    datapath = join(path, makepath(date), prefix + date + "_pci.npy")
    logpath = join(path, makepath(date), prefix + date + "_log.log")
    powpath = join(path, makepath(date), prefix + date + "_pow.npy")
    rfpath = join(path, makepath(date), prefix + date + "_rfi.npy")

    dd = np.load(datapath)
    pd = np.load(powpath)
    powdata.append(np.load(powpath))
    rfdata.append(np.load(rfpath))

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

    pval = []
    for i in range(pd[0,0,:].size):
        pval.append(np.mean(pd[:,:,i]))

    cval.append(np.asanyarray(pval))

    data.append(dd * float(info["Pre-Amp Gain"]))

    xlen = xval.size
    ylen = yval.size
    # plen = pval.size

    inf = np.empty((2,len(info)), dtype='object')
    i = 0
    for k in info:
        inf[0][i] = k
        inf[1][i] = info[k]
        i = i+1

print()
data = np.concatenate((data[0], data[1], data[2]), axis = 2)
rfdata = np.concatenate((rfdata[0], rfdata[1], rfdata[2]), axis = 2)
powdata = np.concatenate((powdata[0], powdata[1], powdata[2]), axis = 2)
pval = np.concatenate((cval[0], cval[1], cval[2]))

savename = join(newpath, newname)
np.savez(savename, d = data, xval = xval, yval = yval, cval = pval, rfi = rfdata, pow = powdata, inf = inf)