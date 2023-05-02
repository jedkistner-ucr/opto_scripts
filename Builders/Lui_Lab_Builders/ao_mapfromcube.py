'''
takes a simple cube and gets a simple map
'''

import numpy as np
from scipy import ndimage
from os.path import join
import sys
from plot_tools import *
from analysis_tools import *
import loader

path = "C:/QMO/Built"
newname = "ao3_deltav"
# name = 'erfu01_sd_bg'
name = 'erfu01_tg_bg'

newname0 = 'Ao1_deltaV'
newname1 = 'Ao1_totalV'

xval, yval, cval, data = loader.load_simplecube(path, name)
#for erfu's maps cval is wavelength and yval is gate (sd is either the other gate or source drain)
#i just want to rearrange these a little
# xval, yval, cval = cval0, xval0, yval0
# data = np.swapaxes(data, 0, 2)
# data = np.swapaxes(data, 2, 1)
# data = np.swapaxes(data, 0, 1)

#delta v and total v
newdata0 = np.zeros((xval.size, cval.size))
newdata1 = np.zeros((xval.size, cval.size))
for i, x in enumerate(xval):
    for j, y in enumerate(yval):
        if y == -x:
            newdata0[i,:] = data[i, j, :]
        if y == x:
            newdata1[i,:] = data[i, j, :]


savename0 = join(path, newname0)
savename1 = join(path, newname1)
np.savez(savename0, d = newdata0, xval = cval, yval = xval)
np.savez(savename1, d = newdata1, xval = cval, yval = xval)