'''
Builds datastack from Erfu's file save system
'''

from os.path import join
from os import rename, listdir
from datetime import date

from plot_tools import *
from analysis_tools import *

path = "E:/Data/Raw/Lui_PL_Data/MW32_SD"

f = listdir(path)

for file in f:
    s = file.split(' ')
    newname = s[0]+"_"+s[-1]
    # print(file)
    # print(newname)
    rename(join(path,file), join(path,newname))