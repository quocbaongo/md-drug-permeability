#!/usr/bin/env python3

import os
import sys
import argparse
import numpy
import math
from uncertainties import ufloat, unumpy
from awh_fep_calc_permeability import loadData

def retrieveHeader(fileName, join=False):

    with open(fileName) as f:
        header = []
        offset = 0

        for line in f:
            line = line.lstrip()
            if len(line) == 0:
                continue
            if line[0] == '@':
                offset += 1
                if 'yaxis' in line:
                    line='@    yaxis label "(kJ mol\\S-1\\N)"\n'
                header += line
            elif line[0] == '#':
                offset += 1
                header += line
            else:
                break

    if join:
        return(''.join(header).split('\n'), offset)
    else:
        return header

def pruneHeader(inHeader, removedColumnIndex, extractColumn = None):

    outHeader = []
    nRemoved = 0
    for line in inHeader.split('\n'):
        #print line
        parts = line.split(' ')
        #print parts
        if len(parts) >= 3 and (len(parts[1]) >= 2 and parts[1][0] == 's') and (len(parts[2]) >= 6 and parts[2] == 'legend'):
            sNumber = int(parts[1][1:])
            if removedColumnIndex - 1 == sNumber - nRemoved or extractColumn != None and (sNumber > 0 and sNumber - nRemoved != extractColumn - 1):
                nRemoved += 1
                continue
            line = '%s s%d %s' % (parts[0], sNumber - nRemoved, ' '.join(parts[2:]))
        line += '\n'
        outHeader.append(line)

    return ''.join(outHeader)

def getAwhTargetError(inFile):

    header, nLines = retrieveHeader(inFile, join=True)
    for line in header:
        if 'AWH metadata: target error' in line:
            parts = line.split()
            return float(parts[9])

def getColumnDimensionality(vals):

    for i in range(1, len(vals)):
        if vals[i] < vals[i-1]:
            return i

    return None

def changeXvgHeaderType(header, newType):

    newHeader = []
    for line in header.split('\n'):
        if '@TYPE' in line:
            line = '@TYPE %s' % newType
        line += '\n'
        newHeader.append(line)

    return ''.join(newHeader)

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

    descriptionStr = 'Makes an -nxy XVG file from a multidimensional file by extracting elements of one specific value in one dimension.'

    parser=argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-i', '--input', type=str,
                        help='Input multidimensional XVG file (XYZ[Z][Z]...).')
    parser.add_argument('-c', '--column', type=int, default=1,
                        help='The column number to extract only specific value from, i.e., in practice remove from the data, e.g. the lambda state column. Default: 1.')
    parser.add_argument('-v', '--data_extraction_index', type=int,
                        help='The entry index to extract from the column specified by the --column option, e.g., 0 if that is a fully interacting lambda state.')
    parser.add_argument('-e', '--extract_column', type=int,
                        help='The column index to keep. If not specified data from all columns (except --column) will be kept, in which case the uncertainty will not be saved in the output.')

    parser.add_argument('--calibration_index', type=int,
                        help='The entry index (in the column to extract) from which to retrieve calibration values, e.g., 30 if that is a fully decoupled lambda state.')
    parser.add_argument('--calibrate_to_value', type=float, nargs='*',
                        help='The target value of --calibration_index. The data values are modified according to the difference between the target value and the actual value. Cannot be used in combination with --calibrate_to_pmf_file. Expects two values (or none), with the second being the uncertainty (1 standard error).')
    parser.add_argument('--calibrate_to_pmf_file', type=str,
                        help='Get the target value of --calibration_index from a PMF xvg file (from a solvation free energy calculation). The xvg file should have the same number of columns (as the extracted data, i.e., without --column) to be able to use the same column number (e.g. bias instead of PMF) as in the extraction and to be able to estimate the uncertainty. Cannot be used in combination with --calibrate_to_value.')

    parser.add_argument('-y', '--symmetrise', action='store_true',
                        help='Symmetrise the output around the origin (X=0) to create a symmetric profile from sampling only one half of the system.')
    parser.add_argument('-r', '--remove_symmetrise_spikes', action='store_true',
                        help='When symmetrising remove spikes from the first and last value (set them to neighbour values).')

    parser.add_argument('-o', '--output', type=str,
                        help='Output file.')

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output
    columnNr = args.column
    valueIndex = args.data_extraction_index
    extractColumn = args.extract_column
    if extractColumn == None:
        extractAll = True
    else:
        extractAll = False
    calibrationIndex = args.calibration_index
    if columnNr == 0:
        xColumn = 1
    else:
        xColumn = 0

    symmetrise = args.symmetrise
    removeSpikes = args.remove_symmetrise_spikes
    if removeSpikes and not symmetrise:
        parser.error('--remove_symmetrise_spikes is only used when --symmetrise is active.')

    calibrateToValue = args.calibrate_to_value
    calibrateToPmfFile = args.calibrate_to_pmf_file

    if calibrateToValue != None and calibrateToPmfFile:
        parser.error('--calibrate_to_value cannot be used in combination with --calibrate_to_pmf_file. Use one or the other.')

    if calibrateToValue and len(calibrateToValue) == 2:
        calibrateToValue = ufloat(calibrateToValue[0], calibrateToValue[1])

    if not inputFile:
        parser.error('An --input file must be specified.')
    if not outputFile:
        parser.error('An --output file must be specified.')
    if valueIndex == None:
        parser.error('A --data_extraction_index must be specified (e.g. the lambda index to keep).')

    if calibrationIndex != None or (calibrateToValue != None or calibrateToPmfFile):
        doCalibration = True
        if calibrationIndex == None or (calibrateToValue == None and calibrateToPmfFile == None):
            parser.error('If calibrating the data --calibration_index and --calibrate_to_value or --calibrate_to_pmf_file must be specified.')
    else:
        doCalibration = False

    data = loadData(inputFile, saveCompressed=False)
    nRows, nCols = data.shape
    x1 = data[:, xColumn]
    x2 = data[:, columnNr]
    nColumn = getColumnDimensionality(x2)
    nX1 = int(nRows / nColumn)
    x1 = numpy.reshape(x1, (nX1, nColumn))
    x2 = numpy.reshape(x2, (nX1, nColumn))
    columnValuesList = []
    for i in range(nCols):
        columnValuesList.append(numpy.reshape(data[:,i], (nX1, nColumn)))

    if extractAll:
        columnData = columnValuesList
    else:
        targetError = getAwhTargetError(inputFile)
        columnData = numpy.empty(nX1, dtype = 'object')
        if doCalibration:
            refIndex = calibrationIndex
        else:
            refIndex = nColumn - valueIndex - 1
        refSampling = columnValuesList[5][:, refIndex]
        valueSampling = columnValuesList[5][:, valueIndex]
        refMetric = columnValuesList[-1][:, refIndex]
        valueMetric = columnValuesList[-1][:, valueIndex]
        invSqrtValueSampling = 1 / numpy.sqrt(valueSampling)
        invSqrtRefSampling = 1 / numpy.sqrt(refSampling)

        ## This needs access to the metric from all input files.
        #uncertainties = numpy.sqrt(targetError*(valueMetric/valueSampling + refMetric/refSampling))
        #for i in range(nX1):
            #columnData[i] = ufloat(columnValuesList[extractColumn][i, valueIndex], uncertainties[i])
        # Assume covariance 1.
        uncertainties = numpy.sqrt(numpy.power(invSqrtValueSampling * targetError, 2) + numpy.power(invSqrtRefSampling * targetError, 2) + 2)
        for i in range(nX1):
            columnData[i] = ufloat(columnValuesList[extractColumn][i, valueIndex], uncertainties[i])

    #print x1
    #print x2
    #print columnValuesList

    if doCalibration and extractColumn != None:
        if calibrateToPmfFile:
            calibrationData = loadData(calibrateToPmfFile, saveCompressed=False)
            calibrateTargetError = getAwhTargetError(calibrateToPmfFile) or 0.25
            valueSampling = calibrationData[0, 4]
            calibrationSampling = calibrationData[-1, 4]
            invSqrtValueSampling = 1 / numpy.sqrt(valueSampling)
            invSqrtCalibrationSampling = 1 / numpy.sqrt(calibrationSampling)
            # Assume covariance 1.
            calibrationUncertainty = numpy.sqrt(numpy.power(invSqrtValueSampling * calibrateTargetError, 2) + numpy.power(invSqrtCalibrationSampling * calibrateTargetError, 2) + 2)
            calibrateToValue = ufloat(calibrationData[-1, extractColumn-1] - calibrationData[0, extractColumn-1], calibrationUncertainty)

        diff = columnValuesList[extractColumn][:,-1] - calibrateToValue
        columnData -= diff

    header = '# Command line:\n'
    header += '# %s\n' % ' '.join(argv)
    inputFileHeader = ''.join(retrieveHeader(inputFile, join=False))

    inputFileHeader = pruneHeader(inputFileHeader, columnNr, extractColumn)

    #print inputFileHeader
    header += inputFileHeader

    outData = (x1[:,valueIndex])
    if extractAll:
        for i in range(nCols):
            if i == xColumn or i == columnNr:
                continue
            outData = numpy.column_stack((outData, columnValuesList[i][:,valueIndex]))
    else:
        outData = numpy.column_stack((outData, unumpy.nominal_values(columnData), unumpy.std_devs(columnData)))
        header = changeXvgHeaderType(header, 'xydy')

    # If there are any NaNs in the output interpolate between existing data points.
    for col in range(outData.shape[1]):
        nans, xNan = nan_helper(outData[:,col])
        outData[nans,col] = numpy.interp(xNan(nans), xNan(~nans), outData[~nans,col])

    if symmetrise:
        nRows, nOutputCols = outData.shape
        if abs(outData[0,0] - 0.0) < 0.0001:
            firstIndex = 1
            symmetriseData = numpy.empty((nRows * 2 - 1, nOutputCols))
        else:
            firstIndex = 0
            symmetriseData = numpy.empty((nRows * 2, nOutputCols))
        nOutputRows = symmetriseData.shape[0]

        symmetriseData[:nRows, 0] = -outData[::-1, 0]
        symmetriseData[nRows:, 0] = outData[firstIndex:, 0]

        if removeSpikes:
            outData[:2,:] = outData[2,:]
            outData[-2:,:] = outData[-3,:]

        for col in range(1, nOutputCols):
            symmetriseData[:nRows, col] = outData[::-1,col]
            symmetriseData[nRows:, col] = outData[firstIndex:, col]

        outData = symmetriseData

    numpy.savetxt(fname=outputFile, X=outData, fmt='%12e', header=header, comments='')

if __name__ == "__main__":
    main(sys.argv)

