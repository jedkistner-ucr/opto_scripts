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
data = data[:,:,10]
# data = data*1e9


save = True
show = True
savename = name + "_signal_lowpass25"
zmin = 0
zmax = 0

#2021_03_25_51
T = 5   #sample time (time per line)
N = 115 #sample count (pixels)
#2022_02_03_12
T = 1   #sample time (time per line)
N = 101 #sample count (pixels)


f0 = 5#6.86 #freq to filter
Q = 2.0 #qaulity

xlen = xval.size
ylen = yval.size
clen = cval.size

axs, figs = makeaxes(3)
ax = []
fig = pl.figure(figsize=figs)
for a in axs:
    ax.append(fig.add_axes(a))

#stack all the data together
# x = np.arange(xlen*ylen)
# val = []
# y = data[0,:]
# for i in range(ylen-1):
#     val.append(data[i+1,:])
# y = np.append(y, val)

#take the fft
# yf = 2*N*ylen *np.abs(fft(y)[:N*ylen//2])
# xf = fftfreq(N*ylen, T*ylen)[:N*ylen//2]
# xf = fftfreq(N*ylen, T*ylen)

#notch filter
b, a = signal.iirnotch(f0, Q, N/T)
b, a = signal.iirfilter(2, f0, fs = N/T, btype='lowpass')
freq, h = signal.freqz(b, a, fs = N/T)
# ax[2].plot(freq * (N/T) / (2*np.pi), np.log(np.abs(h)), c = 'k')
#loop through and calculate the fft of every linetrace
ft = np.zeros((50))
filt_ft = np.zeros((50))
for i in range(ylen):
    x = np.arange(xlen)
    x = xval
    y = data[i,:]
    filtered_y = signal.filtfilt(b, a, y)
    yf = 2*N *np.abs(fft(y)[:N//2])
    ft += yf
    xf = fftfreq(N, T/N)[:N//2]
    filt_yf = 2*N *np.abs(fft(filtered_y)[:N//2])
    filt_ft += filt_yf
    filt_xf = fftfreq(N, T/N)[:N//2]

    ax[0].plot(x, y, c = 'k', alpha = .05)
    ax[2].plot(x, filtered_y, c = 'k', alpha = .05)
    ax[1].plot(xf, yf, c = 'r', alpha = .05)

ft = ft / ylen
filt_ft = filt_ft / ylen
ax[1].plot(xf, yf, c = 'k', alpha = 1)
ax[1].plot(filt_xf, filt_yf, c = 'k', alpha = 1, linestyle = ":")

ax[0].set_xlim(xval[0], xval[-1])
ax[2].set_xlim(xval[0], xval[-1])
ax[1].set_xlim(xf[0], xf[-1])
ax[1].set_ylim(0, np.max(ft)*.1)

ax[0].set_title('Raw curves')
ax[1].set_title('FFT')
ax[2].set_title('Filtered Curves')

ax[2].set_xlabel('Fast Scan Value')
ax[0].set_xlabel('Fast Scan Value')
ax[1].set_xlabel('Frequency (Hz)')

ax[2].set_ylabel('Signal')
ax[0].set_ylabel('Signal')
ax[1].set_ylabel('Amplitude')

if show:
    pl.show()
if save:
    fig.savefig("C:/Jed/Figures/" + savename +".png")