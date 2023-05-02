from builder import build_map
import loader
from analysis_tools import *
from plot_tools import *

#files to load and build
newname = "2023_03_29_sdcube"
d = []
r = []
p = []
cval = []
path = "C:/QMO/Raw Data/2023_03/2023_03_29"
newpath = "C:/QMO/Built"
names  = ["CPS_2023_03_29_97"]
for i in range(20):
    val = 98+i
    names.append("CPS_2023_03_29_" + str(val))

#builds cube taking source drain as the cube axis
for name in names:
    build_map(path, name, newpath)
    xval, yval, data, rf, power, info = loader.load_map(newpath, name, returnpower=True)
    info = make_dict(info)
    d.append(data)
    r.append(rf)
    p.append(power)
    print("Vsd: " + str(info['Source/Drain Start']))
    cval.append(float(info['Source/Drain Start']))
data = np.stack(d, axis = -1)
rf = np.stack(r, axis = -1)
power = np.stack(p, axis = -1)
cval = np.asarray(cval)

savename = join(newpath, newname)
np.savez(savename, d = data, xval = xval, yval = yval, cval = cval, rfi = rf, pow = power, inf = info)