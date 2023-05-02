import numpy as np
from scipy import ndimage
from scipy import interpolate as interp
import scipy.optimize as op
from scipy import fftpack
from scipy.signal import butter, filtfilt
# from curves import diode, exponential, twobody, ln, power, line
from datetime import date
from os.path import join
import sys
sys.path.append("C:/QMO/qmoCode/Toolbox")
from plot_tools import *
from analysis_tools import *
import loader
from builder import build_map
from curves import simplepower



path = "C:/QMO/Built"

name = ["CPS_2022_02_03_26", "CPS_2022_02_03_27", "CPS_2022_02_03_28"]
savename = "CPS_2022_02_03_26_stitch"

xval, yval, data, rf, power, info = [],[],[],[],[],[]
for i in range(len(name)):
    xval1, yval1, data1, rf1, power1, info1 = loader.load_map(path, name[i], returnpower=True)
    xval.append(xval1)
    yval.append(yval1)
    data.append(data1)
    rf.append(rf1)
    power.append(power1)

nx, ny = 1, 1

x = np.linspace(-9, 4, 391)
y = np.linspace(-9, 4, 391)
d = np.zeros((391, 391))
p = np.zeros((391, 391))
r = np.zeros((391, 391))

for i in range(len(data)):
    xstarti = np.searchsorted(x, xval[i][0])
    ystarti = np.searchsorted(y, yval[i][0])

    for n in range(xval[i].size):
        for m in range(yval[i].size):
            if i ==0:
                xindex = xstarti - n
                yindex = ystarti - m
            if i ==1:
                xindex = xstarti - n
                yindex = ystarti - m
            if i ==2:
                xindex = xstarti - n +1
                yindex = ystarti - m +1

            # if d[yindex, xindex] != 0:
            #     print(i)
            #     print(xindex)
            #     print(yindex)
            #     d[yindex, xindex] = 1
            # else:
            d[yindex, xindex] = data[i][m,n]

            p[yindex, xindex] = power[i][m,n]
            r[yindex, xindex] = rf[i][m,n]

# pl.figure()
# pl.imshow(d)
# pl.show()
savename = join(path, savename)
np.savez(savename, d = d, xval = x, yval = y, rfi = r, pow = p, info = info1)

