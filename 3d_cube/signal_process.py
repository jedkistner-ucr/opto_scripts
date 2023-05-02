'''
Takes a built bias hypecube and works with it
'''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy import ndimage, signal
from scipy.fft import fft, fftfreq

from plot_tools import *
from analysis_tools import *
import loader
# from builder import build_map

path = "C:/Jed/Built"
name = "2023_03_29_sdcube" #ao07 sd cube
# name = "2022_02_04_11_sdcube" #d02 sd cube
name = '2022_01_27_0_sdcube'   #sd cube of ao02 at -5.25
# name = '2022_01_27_1_sdcube'  #sd cube of ao02 at -4.48
# name = "2022_02_03_33_powercube" #tg -bg powercube for Ao02
name = "2022_02_03_12_powercube"  #sd_gate powercube ao07
# name = '2021_03_25_51_powercube'
xval, yval, cval, data, rf, powdata, info_ = loader.load_cube(path, name)
# data = data[:,:,25]
# data = data*1e9

T = 1   #sample time (time per line)
N = 101 #sample count (pixels)
f0 = 8#6.86 #freq to filter
Q = 2 #order

newname = name+"_butter_"+str(f0)+"_"+str(Q)
print(newname)

newdata = np.zeros(data.shape)
xlen = xval.size
ylen = yval.size
clen = cval.size

#notch filter
b, a = signal.iirnotch(f0, Q, N/T)
b, a = signal.iirfilter(Q, f0, fs = N/T, btype='lowpass')
freq, h = signal.freqz(b, a, fs = N/T)
# ax[2].plot(freq * (N/T) / (2*np.pi), np.log(np.abs(h)), c = 'k')
#loop through and calculate the fft of every linetrace

for i in range(ylen):
    for q in range(clen):

        y = data[i,:,q]
        filtered_y = signal.filtfilt(b, a, y)
        newdata[i,:,q] = filtered_y


savename = join(path, newname)
np.savez(savename, d = newdata, xval = xval, yval = yval, cval = cval, rfi = rf, pow = powdata, inf = info_)