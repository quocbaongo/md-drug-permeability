#!/usr/bin/env python3

import sys
import numpy
import math
import argparse
from scipy import integrate
from uncertainties import ufloat, unumpy
from awh_fep_calc_permeability import loadData

R = 8.3144598e-3

def calcDiffusionProfile(data):

    # Convert from friction (ps/nm^2) to diffusion (cm^2/s).
    diffusion = numpy.divide(1, (data*100), out=numpy.zeros_like(data), where=data!=0)

    return diffusion

def rolling1DMedian(data, medianWindowSize = 3, minValue = None):

    if medianWindowSize == None:
        medianWindowSize = 3

    nRows, nCols = data.shape
    filteredData = numpy.empty((nRows, nCols))
    filteredData[:] = numpy.nan
    median = numpy.zeros(nRows)

    for r in range(nRows):
        tmpData = data[r,(data[r,:] > minValue)]
        filteredData[r, :len(tmpData)] = tmpData
        #print tmpData
        #print filteredData[r,:]

    #print filteredData

    nElementsOnEachSide = int((medianWindowSize - 1) / 2)

    for r in range(nRows):
        #print data[r-nElementsOnEachSide:r+nElementsOnEachSide, :]
        #print relevantValues[r-nElementsOnEachSide:r+nElementsOnEachSide, :]
        #print filteredData
        median[r] = numpy.nanmedian(filteredData[max(0, r-nElementsOnEachSide):min(nRows, r+nElementsOnEachSide), :])

    #print 'median', median

    return median

def filterData(frictionData, diffFactor = 10, minValue = 6e2, medianWindowSize = 3):

    nRows, nCols = frictionData.shape

    boolMask = numpy.zeros((nRows, nCols), dtype=bool)
    clippedData = numpy.zeros((nRows, nCols))

    rowMedian = rolling1DMedian(abs(frictionData), medianWindowSize, minValue)

    #boolMask = numpy.where(frictionData > rowMax * threshold and frictionData > rowMin, 1, 0)
    for r in range(nRows):
        #boolMask[r,:] = numpy.where(frictionData[r,:] > rowMax[r] * threshold, 1, 0)
        #boolMask[r,:] = numpy.where((frictionData[r,:] > rowMedian[r] / diffFactor) & (frictionData[r,:] < rowMedian[r] * diffFactor), 1, 0)
        boolMask[r,:] = numpy.where(frictionData[r,:] > rowMedian[r] / diffFactor, 1, 0)
        clippedData[r] = numpy.clip(frictionData[r], None, rowMedian[r] * diffFactor)

    return numpy.where(boolMask, clippedData, numpy.nan)

def symmetriseData(data):

    nRows, nCols = data.shape

    # Add a separate column with the mirrored data to count it as separate input.
    symData = numpy.zeros((nRows, nCols * 2))
    symData[:,:nCols] = data[:,:]
    symData[:,nCols:] = data[::-1,:]

    return symData

def dataMeanSemUFloat(data, smoothNPoints = None, useTwoStd = False):

    nRows, nCols = data.shape
    dataMean = numpy.nan_to_num(numpy.nanmean(data, axis=1))

    if smoothNPoints:
        dataMean = smooth(dataMean, smoothNPoints, 'hanning', True)

    dataStd = numpy.nan_to_num(numpy.nanstd(data, axis=1, ddof=1))
    dataCnt = nCols - numpy.count_nonzero(numpy.isnan(data), axis=1)
    dataSem = numpy.divide(dataStd, numpy.sqrt(dataCnt), out=numpy.zeros_like(dataStd), where=dataCnt>1)
    if useTwoStd:
        uncertainty = dataStd * 2
    else:
        uncertainty = dataSem

    out = numpy.empty(nRows, dtype=object)
    for i in range(nRows):
        v = ufloat(dataMean[i], uncertainty[i])
        out[i] = v

    return out

def smooth(x, window_len=11, window='hanning', crop=True):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
        crop: remove points in the beginning and end to keep the length of the input.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.

    Code from: http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')

    if crop:
        half_window = int(window_len/2)
        y = y[half_window:-half_window]

    return y

# From https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]

def main(argv):

    descriptionStr = 'Reads output from friction output from AWH (needs -friction) and calculates the diffusion. The diffusion profile can be smoothened by a Hanning filter.'

    parser=argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-f', '--friction', nargs='+', type=str,
                        help='The friction files from AWH.')
    parser.add_argument('-t', '--temperature', nargs='?', type=float, default=303.15,
                        help='The temperature of the system.')
    parser.add_argument('-c', '--column', type=int,
                        help='Only use this column index (1 is the first friction column).')
    parser.add_argument('-s', '--smooth', nargs='?', type=float, default=None,
                        help='Width (in nm) of Hanning filter to smoothen the diffusion profile.')
    parser.add_argument('-y', '--symmetrise', action='store_true',
                        help='Symmetrise the awh data around the origin (X=0) before processing.')
    parser.add_argument('-r', '--remove_spikes', action='store_true',
                        help='Removes spikes from the first and last value and in the center (set them to neighbour values).')
    parser.add_argument('-o', '--outputbase', type=str,
                        help='The basename of output files to write. Suffixes and extensions will be added.')
    parser.add_argument('--std', action='store_true',
                        help='Report 2 STD instead of 1 SEM as uncertainty.')

    args = parser.parse_args()

    infileNames = args.friction
    if not infileNames:
        parser.error('No AWH friction files specified.')
    outputBase = args.outputbase
    if not outputBase:
        parser.error('Must specify --outputbase.')

    temp = args.temperature
    smoothNm = args.smooth
    symmetrise = args.symmetrise
    removeSpikes = args.remove_spikes
    column = args.column
    useTwoStd = args.std

    if column == 0:
        parser.error('--column cannot be 0. That is the X coordinate.')

    xVals = None

    nFiles = len(infileNames)
    for n, f in enumerate(infileNames):
        data = loadData(f, False)
        if n == 0:
            xVals = data[:,0]
            nRows, nCols = data.shape
            if column != None:
                nDiffusionVals = 1
            else:
                nDiffusionVals = nCols - 1
            frictionData = numpy.zeros((nRows, nFiles * nDiffusionVals))
        if column != None:
            frictionData[:,n] = data[:,column]
        else:
            frictionData[:,n*nDiffusionVals:n*nDiffusionVals+nDiffusionVals] = data[:,1:]

    beta = 1/(R * temp)
    frictionData *= beta

    if removeSpikes:
        #print frictionData[0:4,0], frictionData[nRows/2-1:nRows/2+2,9], frictionData[-4:,0]
        frictionData[:2,:] = frictionData[2,:]
        frictionData[-2:,:] = frictionData[-3,:]
        if nRows % 2 == 1:
            frictionData[nRows//2-1:nRows//2+2,:] = numpy.average(frictionData[[nRows//2-2,nRows//2+3],:], axis=0)
        else:
            frictionData[nRows//2-2:nRows//2+2,:] = numpy.average(frictionData[[nRows//2-3,nRows//2+2],:], axis=0)
        #print frictionData[0:4,0], frictionData[nRows/2-1:nRows/2+2,9], frictionData[-4:,0]

    if smoothNm:
        dataWidthNm = xVals[-1] - xVals[0]
        smoothPoints = int(round(nRows/(dataWidthNm/smoothNm)))
        smoothPoints += (smoothPoints + 1) % 2
    else:
        smoothPoints = None

    # Replace negative values or zeroes with NaN.
    frictionData = numpy.where(frictionData <= 0, numpy.nan, frictionData)

    # If there are any NaNs in the data interpolate between existing data points.
    for col in range(frictionData.shape[1]):
        nans, xNan = nan_helper(frictionData[:,col])
        frictionData[nans,col] = numpy.interp(xNan(nans), xNan(~nans), frictionData[~nans,col])

    if smoothPoints:
        frictionData = filterData(frictionData, diffFactor = 10, medianWindowSize = smoothPoints)
    #print frictionData

    if symmetrise:
        symFrictionData = symmetriseData(frictionData)

    frictionData = dataMeanSemUFloat(frictionData, smoothPoints, useTwoStd)

    if symmetrise:
        symFrictionData = dataMeanSemUFloat(symFrictionData, smoothPoints, useTwoStd)


    diffusionProfile = numpy.empty_like(frictionData)
    diffusionProfile = calcDiffusionProfile(frictionData)


    header = '@TYPE xydy\n'
    header+= '@xaxis label "z (nm)"\r\n'
    header+= '@ legend on\r\n'
    header+= '@ legend box on\r\n'
    header+= '@ legend loctype view\r\n'
    header+= '@ legend 1.1, 0.9\r\n'
    header+= '@ legend length 2\r\n'
    header+= '@ s0 errorbar on\n'

    diffusionFile = outputBase + '.xvg'
    data = numpy.column_stack((xVals, unumpy.nominal_values(diffusionProfile), unumpy.std_devs(diffusionProfile)))
    numpy.savetxt(fname=diffusionFile, X=data, fmt='%12e', header=header, comments='')
    print('Mean diffusion: %s cm^2/s' % numpy.mean(diffusionProfile))

    if symmetrise:
        symDiffusionProfile = numpy.empty_like(symFrictionData)
        symDiffusionProfile = calcDiffusionProfile(symFrictionData)

        symDiffusionFile = outputBase + '_symm.xvg'
        data = numpy.column_stack((xVals, unumpy.nominal_values(symDiffusionProfile), unumpy.std_devs(symDiffusionProfile)))
        numpy.savetxt(fname=symDiffusionFile, X=data, fmt='%12e', header=header, comments='')
        print('Mean diffusion of symmetrised profile: %s cm^2/s' % numpy.mean(symDiffusionProfile))


if __name__ == '__main__':
    main(sys.argv)
