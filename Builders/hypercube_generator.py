'''
stacks cubes into hypercubes. Deoes no alignment or smoothing
'''

import numpy as np
from os.path import join

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Data File Parameters

path = "E:/Data/Raw/2022/2022_01/2022_01_27" #data location where scans are stored
newpath = "C:/Jed/Built"    #new data file save location
newname = "2022_01_27_powermap_noa"                #new data file save name

#The name of the repeating variable in the log file -- this will become zval
# repeatname = "Backgate Start"
repeatname = "Source/Drain Start"

#option to average power data files for cval instead of recorded logfile values
usepower_ascval = True

# Build a list of the data files that you want to allign
#03_25_68 to 03_26_17
# names = ['CPS_2023_03_28_72']
# for i in range(2):
#     num = 73 + i
#     names.append('CPS_2023_03_28_' + str(num))
# for i in range(6):
#     num = 0 + i
#     names.append('CPS_2023_03_29_' + str(num))

# names = ['CPS_2022_01_26_33']
# for i in range(10):
#     num = 34 + i
#     names.append('CPS_2022_01_26_' + str(num))

# names = ['CPS_2022_02_01_49']
# for i in range(5):
#     num = 50 + i
#     names.append('CPS_2022_02_01_' + str(num))
# for i in range(23):
#     num = 0 + i
#     names.append('CPS_2022_02_02_' + str(num))

names = ['CPS_2022_01_27_7']
for i in range(8):
    num = 8 + i
    names.append('CPS_2022_01_27_' + str(num))
for i in range(17):
    num = 0 + i
    names.append('CPS_2022_01_28_' + str(num))
for i in range(7):
    num = 0 + i
    names.append('CPS_2022_01_29_' + str(num))

# Build a list of the data files that you want to allign
#03_24_74 to 03_25_23
# names = ['CPS_2021_03_24_74']
# for i in range(4):
#     num = 75 + i
#     names.append('CPS_2021_03_24_' + str(num))
# for i in range(24):
#     num = 0 + i
#     names.append('CPS_2021_03_25_' + str(num))

print("  ~~  ")
print("Compiling %i data files" % (len(names)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Data Loadup, Axis Generation, Stacking

#compiles images into a single array
#First pull logfile to get xval and yval
logfile = path + "/" + names[0] + "_log.log"
f = open(logfile, "r")
info = {}
info["scan"] = names[0]

for line in f:
        s = line.split(":")
        info[s[0]] = s[1]
f.close()

xval = np.linspace( float(info["Fast Axis Start"]), float(info["Fast Axis End"]), int(info["nx"]))
yval = np.linspace( float(info["Slow Axis Start"]), float(info["Slow Axis End"]), int(info["ny"]))
cval = np.linspace( float(info["Cube Axis Start"]), float(info["Cube Axis End"]), int(info["Cube Scan Number"]))
zval = []

#Zval is pulled from successive log files and rf and data are pulled from successive .npy files
datalist = []
rflist = []
powlist = []
for i in range(len(names)):
    #Pull pci, rfi, and power data files
    d = np.load(join(path, names[i] + "_" + "pci" + ".npy"))
    r = np.load(join(path, names[i] + "_" + "rfi" + ".npy"))
    p = np.load(join(path, names[i] + "_" + "pow" + ".npy"))
    
    datalist.append(d)
    rflist.append(r)
    powlist.append(p)

    logfile = path + "/" + names[i] + "_log.log"
    f = open(logfile, "r")
    for line in f:
        s = line.split(":")
        if s[0] == repeatname:
            zval.append(float(s[1]))
    print(" " + str(len(names)-i))

#Data is stacked
data = np.stack(datalist, axis = -1)
rf = np.stack(rflist, axis = -1)
power = np.stack(powlist, axis = -1)
zval = np.asanyarray(zval)
#data[x, y, c, z]

if usepower_ascval:
    for i in range(cval.size):
        cval[i] = np.mean(power[:,:,i,:])


data = data * float(info["Pre-Amp Gain"])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Saving

savename = join(newpath, newname)
print("saving " + savename)
np.savez(savename, data = data, xval = xval, yval = yval, cval = cval, zval = zval, rfi = rf, pow = power)