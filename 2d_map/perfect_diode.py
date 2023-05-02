'''
Takes a just finished run, builds dataset, creates a map and linecuts
'''

import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
from plot_tools import *
from analysis_tools import *

def diode(x, n, kt, i_0):
    i = i_0*(np.exp(x*n/kt)-1)
    return i
def fermi(x, e, kt, n):
    i = 1 / (np.exp((e-(x*n))/kt)+1)
    return i

xval = np.linspace(-5, 5, 1000)
yval = np.linspace(-5, 5, 1000)
data = np.zeros((yval.size, xval.size))

e = 2
kt = .026 #in eV
n_vsd = 5
n_vg = 1.602*1e-19*1
n_vg = 1

carrier_density = fermi(yval, e, kt, n_vg)

pl.figure()
pl.plot(yval, carrier_density)
pl.show()