import numpy as np
import scipy as sp
from scipy import ndimage
from scipy.optimize import curve_fit
import sys
sys.path.append("C:/Users/Jedki/Documents/analCode/Toolbox")
from plot_tools import make_dict


'''
applies a low pass filter and returns fft and filtered data
'''
def applyGauss(data, filterval):
    newdata = ndimage.gaussian_filter(data, filterval)
    return newdata

def findzeros(data, lowerindex = -1, upperindex = -1):
    zeros = []
    numberTrue = 0
    if lowerindex == -1:
        lowerindex = 0
    if upperindex == -1:
        upperindex = data.size
    for i in range(data.size - 1):
        if i < upperindex and i > lowerindex:
            if data[i] == 0:
                zeros.append(i)
                numberTrue = numberTrue + 1
            elif data[i] < 0 and data[i+1] > 0:
                zeros.append(i)
                numberTrue = numberTrue + 1
            elif data[i] > 0 and data[i+1] < 0:
                zeros.append(i)
                numberTrue = numberTrue + 1
    return zeros, numberTrue

'''
takes a data set and gets the derivative in either x or y
'''
def getdelta(data, dir, degree):
    
    #np.gradient returns the gradient of the input array, for a 2 dimentional array it returns first gradient in rows, then gradient in cols
    
    #choose the direction then iterate over a gradient for each degree
    if dir == "x":
        gradx = data
        for i in range(degree):
            grady, gradx = np.gradient(gradx)
        
        return gradx

    elif dir == "y":
        grady = data
        for i in range(degree):
            grady, gradx = np.gradient(grady)

        return grady

'''
takes xdata and ydata, returns parameters for curve fit
'''
def fit_curves(xdata, ydata, curvetype, startatzero = False, p0 = 0):
    if curvetype == 'power':
        try:
            if p0 == 0:
                par, pcov = curve_fit(power, xdata, ydata)
            else:
                par, pcov = curve_fit(power, xdata, ydata, p0)
            newdata = power(xdata, par[0], par[1])
        except:
            par = [0,0]
            newdata = np.zeros(ydata.shape)
            print('failed')
    elif curvetype == 'diode':
        try:
            if p0 == 0:
                par, pcov = curve_fit(diode, xdata, ydata)
            else:
                par, pcov = curve_fit(diode, xdata, ydata, p0)
            newdata = diode(xdata, par[0], par[1])
        except:
            par = [0,0]
            newdata = np.zeros(ydata.shape)
            print('failed')
    elif curvetype == 'exponential':
        try:
            if p0 == 0:
                par, pcov = curve_fit(exponential, xdata, ydata)
            else:
                par, pcov = curve_fit(exponential, xdata, ydata, p0)
            newdata = exponential(xdata, par[0], par[1], par[2])
        except:
            par = [0,0,0]
            newdata = np.zeros(ydata.shape)
            print('failed')
    elif curvetype == 'sigmoid':
        try:
            if p0 == 0:
                par, pcov = curve_fit(sigmoid, xdata, ydata)
            else:
                par, pcov = curve_fit(sigmoid, xdata, ydata, p0)
            newdata = sigmoid(xdata, par[0], par[1], par[3])
        except:
            par = [0,0]
            newdata = np.zeros(ydata.shape)
            print('failed')
    elif curvetype == 'neggauss':
        try:
            if p0 == 0:
                par, pcov = curve_fit(neggauss, xdata, ydata)
            else:
                par, pcov = curve_fit(neggauss, xdata, ydata, p0)
            newdata = neggauss(xdata, par[0], par[1], par[2])
        except:
            par = [0,0,0]
            newdata = np.zeros(ydata.shape)
            # print('failed')

    elif curvetype == 'negloren':
        try:
            if p0 == 0:
                par, pcov = curve_fit(negloren, xdata, ydata)
            else:
                par, pcov = curve_fit(negloren, xdata, ydata, p0)
            newdata = negloren(xdata, par[0], par[1], par[2])
        except:
            par = [0,0,0]
            newdata = np.zeros(ydata.shape)
            # print('failed')
    elif curvetype == 'linear':
        try:
            if p0 == 0:
                par, pcov = curve_fit(line, xdata, ydata)
            else:
                par, pcov = curve_fit(line, xdata, ydata, p0)
            newdata = line(xdata, par[0], par[1])
        except:
            par = [0,0]
            newdata = np.zeros(ydata.shape)
            # print('failed')

    elif curvetype == 'twobody':
        if startatzero:
            ydata[0] = 0
        try:
            if p0 == 0:
                par, pcov = curve_fit(twobody, xdata, ydata)
            else:
                par, pcov = curve_fit(twobody, xdata, ydata, p0)
            newdata = twobody(xdata, par[0], par[1], par[2])
        except:
            par = [0,0,0]
            newdata = np.zeros(ydata.shape)

    elif curvetype == 'twobody_line':
        try:
            if p0 == 0:
                par, pcov = curve_fit(twobody_line, xdata, ydata)
            else:
                par, pcov = curve_fit(twobody_line, xdata, ydata, p0)
            newdata = twobody_line(xdata, par[0], par[1], par[2], par[3])
        except:
            par = [0,0,0, 0]
            newdata = np.zeros(ydata.shape)
            # print('failed')
    # if par[0] < 0:
    #     par = [0, 0, 0]
    #     newdata = np.zeros(ydata.shape)
    return newdata, par
# def fitpoly(xdata, ydata, dim):
#     try:
#         par, res = sp.polyfit(xdata, ydata, dim)
#         x = 0
#         for i in range(dim + 1):
#             x = x + par[i] ** (i+1)
#     except:
#         par = []
#         for i in range(dim+1):
#             par.append(0)
#         newdata = data
    
#     return newdata, par

'''
Takes a 2d set of data and an interval and retuns data cuts and values of cuts
dir 0 = x dir 1 = y
'''
def get_slices(data, xrange, yrange, step = 1, dir = 0):
    slices = []
    zrange = []
    if dir == 0:
        for i in range(0, yrange.size, step):
            slices.append(data[ i , : ])
            zrange.append(yrange[i])
    elif dir == 1:
        for i in range(0, xrange.size, step):
            slices.append(data[ :, i ])
            zrange.append(xrange[i])

    return slices, zrange
'''
Takes two supersets and makes them the same size by getting rid of non mathing data
Right now this is built such that superset1 should be larger than superset2
'''
def trimtemptocompare(x1, y1, d1, r1, p1, i1, x2, y2, d2, r2, p2, i2):
    margin = 2      #temperature margin to keep dataset
    indexdel = []   #indexes that will be deleted
    
    in1 = np.empty(i1[:,0,0].shape, dtype=object)
    in2 = np.empty(i2[:,0,0].shape, dtype=object)
    for i in range(i1[:,0,0].size):
        in1[i] = make_dict(i1[i,:,:])
    for i in range(i2[:,0,0].size):
        in2[i] = make_dict(i2[i,:,:])
    

    #Looping through larger data
    for i in range(x1[:,0].size):
        temp = in1[i]['temp']
        exists = False
        #Check the smaller superset to see if the same temperature exists there
        for q in range(x2[:,0].size):
            temp2 = in2[q]['temp']
            if temp2 < (temp + margin) and temp2 > (temp - margin):
                exists = True
        #If it doesnt exists add the index to the delete list
        if exists == False:
            indexdel.append(i)
        
    #Delete those entries in reverse order
    for i in range(len(indexdel)):
        n = indexdel[ len(indexdel) -1 - i]
        x1 = np.delete(x1, n, axis = 0)
        y1 = np.delete(y1, n, axis = 0)
        d1 = np.delete(d1, n, axis = 0)
        r1 = np.delete(r1, n, axis = 0)
        p1 = np.delete(p1, n, axis = 0)
        i1 = np.delete(i1, n, axis = 0)
    
    return x1, y1, d1, r1, p1, i1
        


    