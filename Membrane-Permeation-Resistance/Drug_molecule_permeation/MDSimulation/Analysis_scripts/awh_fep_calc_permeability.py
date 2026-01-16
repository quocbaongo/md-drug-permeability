#!/usr/bin/env python3

import argparse
import math
import numpy
import os
from scipy import integrate, signal, stats
from sklearn.cluster import MeanShift, estimate_bandwidth
from statsmodels.stats.weightstats import DescrStatsW
import sys
from uncertainties import ufloat, unumpy, umath, correlation_matrix, correlated_values_norm

R = 8.3144598e-3

PROGRESS_BAR_LENGTH = 20

def loadData(fname, saveCompressed=True):

    loadedNpData = False

    #print fname

    if fname[-4:] == '.npz':
        data = numpy.load(fname)['data']
        loadedNpData = True
    elif fname[-4:] == '.npy':
        data = numpy.load(fname)
        loadedNpData = True
    else:
        if os.path.isfile(fname+'.npz'):
            data = numpy.load(fname+'.npz')['data']
            loadedNpData = True
        elif os.path.isfile(fname+'.npy'):
            data = numpy.load(fname+'.npy')
            loadedNpData = True
        else:
            try:
                (xHeader, xOffset) = retrieveHeader(fname)
                data = numpy.loadtxt(fname, skiprows=xOffset, comments='@')
            except Exception:
                print('Error opening file: %s' % fname)
                raise

    # There can be duplicate rows when appending to a pull file. Remove them.
    tmp = data.ravel().view(numpy.dtype((numpy.void, data.dtype.itemsize*data.shape[1])))
    _, unique_idx = numpy.unique(tmp, return_index=True)
    data = data[numpy.sort(unique_idx)]

    if not loadedNpData and saveCompressed:
        numpy.savez_compressed(fname+'.npz', data=data)

    return data

def retrieveHeader(fileName):

    with open(fileName) as f:
        header = []
        offset = 0

        for line in f:
            line = line.strip()
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

    return(''.join(header).split('\n'), offset)


## From https://stackoverflow.com/questions/44617671/how-to-propagate-uncertainty-to-a-user-defined-function-using-python
## Results are the same as numpy.trapz with ufloats.
#def utrapz(y, x):

    #deltaX = x[1:]-x[:-1]
    #avgY = (y[1:]+y[:-1])/2
    #return numpy.sum(deltaX*avgY)

def calcPermeationResistance(xData, pmfData, diffusionData, beta):

    nX = xData.shape[0]

    if isinstance(pmfData[0], float):
        resistProf = numpy.empty(nX)
        resistProf = numpy.exp(pmfData*beta)
    else:
        resistProf = numpy.empty(nX, dtype=object)
        for i in range(nX):
            resistProf[i] = umath.exp(pmfData[i]*beta)

    resistProf /= diffusionData

    # Convert xData from nm to cm
    permeationResistance = numpy.trapz(resistProf, xData * 1e-7)

    return permeationResistance, resistProf

# Calculated according to eq 7 in Das et al. Soft Matter, 2009, 5, 4549-4555.
def calcCrossingTimeOld(xData, pmfData, diffusionData, beta, nLayers=1, diffusionLayerThickness=0, diffusionLayerDiffusionCoefficient=2.3e-5, diffusionLayerDGOffset=0.0):

    if nLayers > 1:
        xMin = xData[0]
        xMax = xData[-1]
        nX = xData.shape[0]
        xDataLocal = numpy.linspace(xMin * nLayers, xMax * nLayers, nX * nLayers)
        pmfDataLocal = numpy.repeat(pmfData, nLayers)
        diffusionDataLocal = numpy.repeat(diffusionData, nLayers)

    else:
        xDataLocal = xData
        pmfDataLocal = pmfData
        diffusionDataLocal = diffusionData

    if diffusionLayerThickness > 0:
        xDataStep = xDataLocal[-1] - xDataLocal[-2]
        xDataLocal = numpy.append(xDataLocal, [xDataLocal[-1]+xDataStep, xDataLocal[-1]+xDataStep+diffusionLayerThickness])
        diffusionDataLocal = numpy.append(diffusionDataLocal, [diffusionLayerDiffusionCoefficient, diffusionLayerDiffusionCoefficient])
        pmfDataLocal = numpy.append(pmfDataLocal, [diffusionLayerDGOffset, diffusionLayerDGOffset])

        #print xDataLocal[-10:]
        #print diffusionDataLocal[-10:]
        #print pmfDataLocal[-10:]

    nX = xDataLocal.shape[0]

    xDataCm = xDataLocal*1e-7

    if isinstance(pmfDataLocal[0], float):
        crossingTime = numpy.empty(nX)
        firstTermProfile = numpy.exp(-pmfDataLocal * beta)
        expr = numpy.exp(pmfDataLocal * beta) / diffusionDataLocal
    else:
        firstTermProfile = numpy.empty(nX, dtype=object)
        permeationResistanceProfile = numpy.empty(nX, dtype=object)
        expr = numpy.empty(nX, dtype=object)
        for i in range(nX):
            firstTermProfile[i] = umath.exp(-pmfDataLocal[i] * beta) # Unitless
            permeationResistanceProfile[i] = umath.exp(pmfDataLocal[i] * beta) # Unitless (will become s/cm^2)
        permeationResistanceProfile /= diffusionDataLocal # s/cm^2

    #print 'permeationResistance', permeationResistanceProfile[-10:]
    # The integration will have to be reversed. integratedPermeationResistanceProfile[0] should be the integral of permeationResistanceProfile[0:],
    # integratedSecondTermProfile[10] should be the integral of permeationResistanceProfile[10:]
    integratedPermeationResistanceProfile = integrate.cumtrapz(permeationResistanceProfile[::-1], -xDataCm[::-1], initial=0)[::-1]
    #print 'integrated permeationResistance', integratedPermeationResistanceProfile[:10], integratedPermeationResistanceProfile[-10:]

    for i in range(nX):
        expr[i] = firstTermProfile[i] * integratedPermeationResistanceProfile[i]
        #print i, expr[i], firstTermProfile[i], integratedPermeationResistanceProfile[i]

    #print 'x (cm)', xDataCm[-10:]
    crossingTime = numpy.trapz(expr, xDataCm)
    #print 'cumTime', integrate.cumtrapz(expr, xDataCm)

    return crossingTime

# Calculated according to eq. 21 in Votapka et al. J. Chem. Phys. B, 2016, 120, 8606-8616.
def calcCrossingTime(pmfXVals, pmfData, permeabilityCoefficient, beta, nLayers=1, diffusionLayerThickness=0, diffusionLayerDiffusionCoefficient=2.3e-5, diffusionLayerDGOffset=0.0):

    nX = pmfXVals.shape[0]

    if isinstance(pmfData[0], float):
        prof = numpy.empty(nX)
        prof = numpy.exp(-pmfData*beta)
    else:
        prof = numpy.empty(nX, dtype=object)
        for i in range(nX):
            prof[i] = umath.exp(-pmfData[i]*beta)

    xDataCm = pmfXVals*1e-7

    permeationResistanceCoefficient = 1/permeabilityCoefficient
    if diffusionLayerThickness > 0:
        permeationResistanceCoefficient += diffusionLayerThickness*1e-7 * (math.exp(diffusionLayerDGOffset*beta)/diffusionLayerDiffusionCoefficient)

    crossingTime = (numpy.trapz(prof, xDataCm) * nLayers + diffusionLayerThickness*1e-7 * diffusionLayerDGOffset) / (2 * (1/permeationResistanceCoefficient) / nLayers)
    crossingTime *= 3600

    return crossingTime

def calculatePmfReferenceValue(pmf, beta):

    n = pmf.shape[0]
    unitlessPmf = numpy.empty(n, dtype=object)
    for i in range(n):
        unitlessPmf[i] = umath.exp(-pmf[i]*beta)

    dimensionlessRefValue = numpy.mean(unitlessPmf)
    refValue = -unumpy.log(dimensionlessRefValue)/beta

    return refValue

def arrayToUfloatArray(array, errorSrc = None):

    nRows, nCols = array.shape

    ufloatArray = numpy.empty(nRows, dtype=object)

    for i in range(nRows):
        val = ufloat(array[i,0], array[i,1], errorSrc)
        ufloatArray[i] = val

    return ufloatArray

def errorSourcesToStr(sum_value):
    outStr = ''
    tags = set(k.tag for k in list(sum_value.error_components().keys()))
    for tag in tags:
        if tag != None:
            tot_error = math.sqrt(sum(err**2 for (var, err) in list(sum_value.error_components().items()) if var.tag == tag))
            outStr += '%s: %s\n' % (tag, tot_error)

    return outStr

def getFullCoordinateRange(dataArrayList):

    nVals = numpy.max([a.shape[0] for a in dataArrayList])
    minVal = numpy.min([a[:,0] for a in dataArrayList])
    maxVal = numpy.max([a[:,0] for a in dataArrayList])

    return numpy.linspace(minVal, maxVal, nVals)

def interpolateToCoordinates(dataArrayList, coords):

    dataArray = numpy.empty((len(coords), len(dataArrayList) + 1))
    dataArray[:, 0] = coords
    for i, data in enumerate(dataArrayList):
        dataArray[:, i+1] = numpy.interp(coords, data[:, 0], data[:, 1])

    return dataArray

def matchDataToCoordinates(dataArrayList):

    coords = getFullCoordinateRange(dataArrayList)

    dataArray = interpolateToCoordinates(dataArrayList, coords)

    return dataArray

def getPeakClusterPositions(peaks, maxDistanceForMatching=0.1):

    zippedPeaks = numpy.array(list(zip(peaks,numpy.zeros(len(peaks)))))
    ms = MeanShift(bandwidth=maxDistanceForMatching/2.0, bin_seeding=True, cluster_all=True, min_bin_freq=1)
    ms.fit(zippedPeaks)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = numpy.unique(labels)
    n_clusters = len(labels_unique)

    # For some reason cluster_centers are not always in the center.
    peakCenters = numpy.empty(n_clusters)
    for k in range(n_clusters):
        my_members = labels == k
        peakCenters[k] = numpy.mean(my_members)

    return peakCenters


def fitToPeaksInData(dataArray, maxDistanceForMatching=0.1, minimalDistanceBetweenPeaks=1.0, heightMeanDifference=8.0):

    peaks = []
    numPeaksPerNm = 1.0/(dataArray[1, 0] - dataArray[0, 0])
    for data in dataArray[:, 1]:
        dataMean = numpy.mean(data)
        maxima=(signal.find_peaks(data, height=dataMean+heightMeanDifference, distance=numPeaksPerNm * minimalDistanceBetweenPeaks))
        minima=(signal.find_peaks(-data, height=-dataMean+heightMeanDifference, distance=numPeaksPerNm * minimalDistanceBetweenPeaks))
        peaks += maxima.tolist() + minima.tolist()

    peaks = numpy.array(peaks)
    #print peaks
    peakAveragePositions = getPeakClusterPositions(peaks, maxDistanceForMatching)

def averagePmfs(dataArray, pmfWeights=None, useTwoStd=False):

    stats = DescrStatsW(numpy.transpose(dataArray[:,1:]), weights=pmfWeights, ddof=1)

    nRows = dataArray.shape[0]
    xAndResults = numpy.empty((nRows, 3), dtype=object)
    xAndResults[:,0] = dataArray[:,0]
    xAndResults[:,1] = stats.mean
    if useTwoStd:
        xAndResults[:,2] = stats.std * 2
    else:
        xAndResults[:,2] = stats.std_mean

    return xAndResults

def main(argv):

    descriptionStr = 'Calculates the permeability from a PMF and a local diffusion coefficient profile.'

    parser=argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-p', '--pmf', type=str, nargs='+',
                        help='XVG file(s) containing the pmf. If multiple files are submitted, the uncertainty will be the estimated standard error from bootstrapping.')
    parser.add_argument('-w', '--pmfweights', type=float, nargs='+',
                        help='If multiple pmf XVG files are used as input, they can be weighted by factors specified by this option. If any --pmfweights are specified they must be in the same order, and the same number, as the files specified by --pmf.')
    parser.add_argument('-d', '--diffusion', type=str,
                        help='XVG file containing the diffusion profile. This cannot be multiple files, but should already be calculated from multiple input friction files using, e.g., awh_diffusion.py.')
    parser.add_argument('-t', '--temperature', nargs='?', type=float, default=303.15,
                        help='The temperature of the system. Default 303.15 K')
    parser.add_argument('-f', '--diffusion_factor', type=float, default=1.0,
                        help='Factor to multiply the diffusion, e.g., if it is known to be too high or too low.')
    parser.add_argument('-z', '--auto_zero_up', action='store_true',
                        help='Automatically shift the PMF up so that its minimum (including the reference solvent) is zero.')
    parser.add_argument('-zz', '--auto_zero_full', action='store_true',
                        help='Automatically shift the PMF so that its minimum is zero.')
    parser.add_argument('--auto_zero_full_after_first_layer', action='store_true',
                        help='Automatically shift the PMF so that its minimum is zero for all layers except the first.')
    parser.add_argument('-zzz', '--auto_ref', action='store_true',
                        help='Automatically shift the PMF related to the reference that is the mean of e^(-beta*deltaG).')
    parser.add_argument('-l', '--layers', nargs='+', type=float, default=[1],
                        help='Calculate the permeability coefficient (and lag time) using this number of layers of the system. This can either be a single integer or a pair of floats, with the second being the SEM.')
    #parser.add_argument('--diffusion_layer_thickness', type=float, default=0,
                        #help='Add an additional layer with pure diffusion of this thickness (nm).')
    #parser.add_argument('--diffusion_layer_diffusion_coefficient', type=float, default=2.3e-5,
                        #help='The diffusion coefficient of the molecule in the diffusion layer (cm2/s).')
    #parser.add_argument('--diffusion_layer_dg_offset', type=float,
                        #help='The difference in free energy between the diffusion layer and the solvent (default 0.0). Cannot be used together with --diffusion_layer_dg_awh_file.')
    #parser.add_argument('--diffusion_layer_dg_awh_file', type=str,
                        #help='Results from AWH decoupling in a medium similar to the diffusion layer (usually water). Cannot be used together with --diffusion_layer_dg_offset.')
    parser.add_argument('--std', action='store_true',
                        help='Report 2 STD instead of 1 SEM as uncertainty.')
    parser.add_argument('-o', '--outputbase', type=str,
                        help='The base name of the files to write results (txt), PMF (xvg), diffusion profile (xvg) and permeation resistance (xvg) output.')

    args = parser.parse_args()

    pmfFiles = args.pmf
    diffusionFile = args.diffusion
    outputBase = args.outputbase
    temp = args.temperature
    beta = 1/(R*temp)

    if not pmfFiles or not diffusionFile:
        parser.error('Must specify --pmf and --diffusion.')

    diffusionFactor = args.diffusion_factor
    autoZeroUp = args.auto_zero_up
    autoZeroFull = args.auto_zero_full
    autoZeroAfterFirstLayer = args.auto_zero_full_after_first_layer
    autoRef = args.auto_ref
    refereceOptionCount = 1 if autoZeroUp else 0 + 1 if autoZeroFull else 0 + 1 if autoZeroAfterFirstLayer else 0 + 1 if autoRef else 0
    if refereceOptionCount > 1:
        parser.error('Cannot combine multiple referencing options.')
    pmfWeights = args.pmfweights
    if pmfWeights:
        if len(pmfWeights) != len(pmfFiles):
            parser.error('When specifying weights for the input PMF files the weights must be as many as the input files.')
    else:
        pmfWeights = None

    nLayers = args.layers
    if len(nLayers) == 1:
        nLayers = int(nLayers[0])
    elif len(nLayers) == 2:
        nLayers = ufloat(nLayers[0], nLayers[1])
    else:
        parser.error('nLayers can either be a single integer or a pair of floats where the second one is the SEM of the first one.')

    # No zeroing will be done if there is not more than one layer.
    if autoZeroAfterFirstLayer and not nLayers > 1:
        autoZeroAfterFirstLayer = False

    # Assume 10% uncertainty in epidermis/full skin thickness and diffusion coefficient therein.
    #diffusionLayerThickness = args.diffusion_layer_thickness
    #if diffusionLayerThickness:
        #diffusionLayerThickness = ufloat(diffusionLayerThickness, diffusionLayerThickness*0.1)
    #diffusionLayerDiffusionCoefficient = args.diffusion_layer_diffusion_coefficient
    #if diffusionLayerDiffusionCoefficient:
        #diffusionLayerDiffusionCoefficient = ufloat(diffusionLayerDiffusionCoefficient, diffusionLayerDiffusionCoefficient*0.1)
    #diffusionLayerDGOffset = args.diffusion_layer_dg_offset
    #diffusionLayerAwhFile = args.diffusion_layer_dg_awh_file
    useTwoStd = args.std

    ## FIXME: diffusionLayerAwhFile is not enabled yet.
    #if diffusionLayerAwhFile != None:
        #parser.error('diffusion_layer_dg_awh_file is not enabled yet.')

    #if diffusionLayerDGOffset != None and diffusionLayerAwhFile != None:
        #parser.error('--diffusion_layer_dg_offset cannot be used together with --diffusion_layer_dg_awh_file as both represent the same thing.')
    #if diffusionLayerDGOffset == None and diffusionLayerAwhFile == None:
        #diffusionLayerDGOffset = 0.0

    if len(pmfFiles) == 1:
        dataArray = loadData(pmfFiles[0], False)
        #dataHeader, nHeaderLines = retrieveHeader(pmfFile)

    else:
        dataArrayList = []
        for pmfFile in pmfFiles:
            dataArrayList.append(loadData(pmfFile, False))

        dataArray = matchDataToCoordinates(dataArrayList)

        #dataArray = fitToPeaksInData(dataArray, 0.1, 1.0, 8.0)
        #dataArray = averageAndUncertaintyFromBootstrap(dataArray)

        #dataArray = loadData(pmfFile, False)

        dataArray = averagePmfs(dataArray, pmfWeights, useTwoStd)

    diffusionArray = loadData(diffusionFile, False)

    pmfXVals = dataArray[:,0]
    diffusionXVals = diffusionArray[:, 0]

    #print dataArray

    nRows = pmfXVals.shape[0]
    #print nRows, math.pow(nRows, 1.0/4)

    if dataArray.shape[1] > 2:
        # Compensate for artificially low uncertainty due to correlated errors
        dataArray[:, 2] *= math.pow(2*nRows, 1.0/4)
    else:
        # Otherwise add a column and set an arbitrarily high error.
        dataArray = numpy.c_[dataArray, numpy.ones(nRows)]
        dataArray[:, 2] *= 10
    pmf = arrayToUfloatArray(dataArray[:,1:], 'PMF')
    pmfMin = numpy.amin(pmf)
    if autoZeroUp or autoZeroFull or autoZeroAfterFirstLayer:
        uncertainty = pmfMin.std_dev
        if autoZeroUp:
            # Probability of the lowest value being below 0.
            probabilityFactor = stats.norm.cdf(0, pmfMin.nominal_value, uncertainty)
            shiftValue = min(0, pmfMin.nominal_value)
        else:
            # Probability of the lowest value not being 0.
            probabilityFactor = max(0, 1-stats.norm.pdf(0, pmfMin.nominal_value, uncertainty))
            shiftValue = pmfMin.nominal_value
        if autoZeroAfterFirstLayer:
            pmfOrig = pmf.copy()
        uncertainty *= probabilityFactor
        pmf -= ufloat(shiftValue, uncertainty / math.pow(2*nRows, 1.0/4))
        #diffusionLayerDGOffset -= pmfMin.nominal_value # TODO: Check why the diffusionLayer dG must be 0 to give reasonable results.
    elif autoRef:
        pmfMin = calculatePmfReferenceValue(pmf, beta)
        print('Reference value: %s kJ/mol' % pmfMin)
        pmf -= ufloat(pmfMin.nominal_value, pmfMin.std_dev / math.pow(2*nRows, 1.0/4))

    diffusionProfile = arrayToUfloatArray(diffusionArray[:,1:] * diffusionFactor, 'Diffusion')

    diffusionNominalVals = numpy.interp(unumpy.nominal_values(pmfXVals), unumpy.nominal_values(diffusionXVals), unumpy.nominal_values(diffusionProfile))
    diffusionStdevs = numpy.interp(unumpy.nominal_values(pmfXVals), unumpy.nominal_values(diffusionXVals), unumpy.std_devs(diffusionProfile))

    diffusionProfile = numpy.empty(nRows, dtype='object')

    for i in range(nRows):
        diffusionProfile[i] = ufloat(diffusionNominalVals[i], diffusionStdevs[i], 'Diffusion')

    permeationResistance, permeationResistanceProfile = calcPermeationResistance(pmfXVals, pmf, diffusionProfile, beta)
    if autoZeroAfterFirstLayer:
        permeationResistanceOrig, permeationResistanceProfileOrig = calcPermeationResistance(pmfXVals, pmfOrig, diffusionProfile, beta)

    # Also save a PMF and a diffusion plot with error bars (xydy cannot be in the same xvg file as xy).
    header = '@TYPE xydy\n'
    header+= '@xaxis label "z (nm)"\r\n'
    header+= '@ legend on\r\n'
    header+= '@ legend box on\r\n'
    header+= '@ legend loctype view\r\n'
    header+= '@ legend 1.1, 0.9\r\n'
    header+= '@ legend length 2\r\n'
    header+= '@ s0 errorbar on\n'

    if outputBase:
        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(pmf), unumpy.std_devs(pmf) / math.pow(2*nRows, 1.0/4)))
        numpy.savetxt(fname='%s_pmf.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(diffusionProfile), unumpy.std_devs(diffusionProfile)))
        numpy.savetxt(fname='%s_diffusion.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(permeationResistanceProfile), unumpy.std_devs(permeationResistanceProfile)))
        numpy.savetxt(fname='%s_perm_resistance.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

    ## Include uncertainties in the minimum
    #if (autoZeroUp and pmfMin < 0) or autoZeroFull:
        #pmf += pmfMin.nominal_value
        #pmf -= pmfMin

    permeationResistanceUnc = permeationResistance.std_dev
    permeationResistance = ufloat(permeationResistance.nominal_value, permeationResistanceUnc)
    p = 3600/permeationResistance
    #crossingTime = calcCrossingTimeOld(pmfXVals, pmf, diffusionProfile, beta)
    #print 'old Crossing Time', crossingTime
    crossingTime = calcCrossingTime(pmfXVals, pmf, p, beta)
    if autoZeroAfterFirstLayer:
        permeationResistanceUncOrig = permeationResistanceOrig.std_dev
        permeationResistanceOrig = ufloat(permeationResistanceOrig.nominal_value, permeationResistanceUncOrig)
        permeationResistanceFull = permeationResistanceOrig + (permeationResistance * (nLayers - 1))
        pOrig = 3600/permeationResistanceOrig
        pFull = 3600/permeationResistanceFull
        layersCrossingTime = calcCrossingTime(pmfXVals, pmf, p, beta, nLayers-1)
        crossingTimeOrig = calcCrossingTime(pmfXVals, pmfOrig, pOrig, beta)
        crossingTimeFull = layersCrossingTime + crossingTimeOrig

    #if diffusionLayerThickness > 0:
        #diffusionPermeationResistance = diffusionLayerThickness * 1e-7 * (math.exp(diffusionLayerDGOffset*beta)/diffusionLayerDiffusionCoefficient)
        ##totalCrossingTime = calcCrossingTimeOld(pmfXVals, pmf, diffusionProfile, beta, nLayers, diffusionLayerThickness, diffusionLayerDiffusionCoefficient, diffusionLayerDGOffset)

    if autoZeroAfterFirstLayer:
        print('Using %s layers of barrier (zeroing the PMF after the first layer):' % nLayers)
        print('  PermeationResistance: %s s/cm (%s h/cm)' % (permeationResistanceFull, permeationResistanceFull / 3600))
        print('  P: %s cm/h' % pFull)
        print('  logP: %s' % umath.log10(pFull))
        print('  Crossing time: %s s' % crossingTimeFull)
    else:
        print('PermeationResistance: %s s/cm (%s h/cm)' % (permeationResistance, permeationResistance / 3600))
        print('P: %s cm/h' % p)
        #print 'P error sources:'
        #print errorSourcesToStr(p),
        #print ''
        print('logP: %s' % umath.log10(p))
        print('Crossing time: %s s' % crossingTime)
        if nLayers > 1:
            print('Using %s layers of barrier:' % nLayers)
            print('  PermeationResistance: %s s/cm (%s h/cm)' % (permeationResistance * nLayers, permeationResistance * nLayers / 3600))
            print('  P: %s cm/h' % (p/nLayers))
            print('  logP: %s' % umath.log10(p/nLayers))
            layersCrossingTime = calcCrossingTime(pmfXVals, pmf, p, beta, nLayers)
    #        print '  Crossing time:', (crossingTime * nLayers * nLayers), 's'
            print('  Crossing time: %s s' % (layersCrossingTime))
        #if diffusionLayerThickness > 0:
            #print('Using %s layers of barrier and %s nm of pure diffusion with diffusion coefficient %s cm2/s:' % (nLayers, diffusionLayerThickness, diffusionLayerDiffusionCoefficient))
            #print('  P:', (3600/(permeationResistance*nLayers+diffusionPermeationResistance)), 'cm/h')
            #print('  logP:', umath.log10((3600/(permeationResistance*nLayers+diffusionPermeationResistance))))
            #totalCrossingTime = calcCrossingTime(pmfXVals, pmf, p, beta, nLayers, diffusionLayerThickness, diffusionLayerDiffusionCoefficient, diffusionLayerDGOffset)
            #print('  Crossing time: %s s' % totalCrossingTime)

    if outputBase:
        with open(outputBase+'.txt', 'w') as f:
            f.write(' '.join(argv) + '\n\n')
            if autoZeroAfterFirstLayer:
                f.write('Using %s layers of barrier (zeroing the PMF after the first layer):' % nLayers)
                f.write('PermeationResistance: %s s/cm (%s h/cm)\n' % (permeationResistanceFull, permeationResistanceFull / 3600))
                f.write('P: %s cm/h\n' % pFull)
                f.write('logP: %s\n' % umath.log10(pFull))
                f.write('Crossing time: %s s\n' % crossingTimeFull)
            else:
                if autoRef:
                    f.write('Reference value: %s kJ/mol\n' % pmfMin)
                f.write('PermeationResistance: %s s/cm (%s h/cm)\n' % (permeationResistance, permeationResistance / 3600))
                f.write('P: %s cm/h\n' % p)
                #f.write('P error sources:\n')
                #srcStr = errorSourcesToStr(p)
                #f.write(srcStr)
                #f.write('\n')
                f.write('logP: %s\n' % umath.log10(p))
                f.write('Crossing time: %s s\n' % crossingTime)
                if nLayers > 1:
                    f.write('Using %s layers of barrier:\n' % nLayers)
                    f.write('  PermeationResistance: %s s/cm (%s h/cm)\n' % (permeationResistance * nLayers, permeationResistance * nLayers / 3600))
                    f.write('  P: %s cm/h\n' % (p/nLayers))
                    f.write('  logP: %s\n' % umath.log10(p/nLayers))
                    f.write('  Crossing time: %s s\n' % (layersCrossingTime))
                #if diffusionLayerThickness > 0:
                    #f.write('Using %s layers of barrier and %s nm of pure diffusion with diffusion coefficient %s cm2/s:\n' % (nLayers, diffusionLayerThickness, diffusionLayerDiffusionCoefficient))
                    #f.write('  P: %s cm/h\n' % (3600/(permeationResistance*nLayers+diffusionPermeationResistance)))
                    #f.write('  logP: %s\n' % umath.log10((3600/(permeationResistance*nLayers+diffusionPermeationResistance))))
                    ##f.write('  Crossing time: %s s\n' % ((crossingTime * nLayers * nLayers) + (math.pow(diffusionLayerThickness * 1e-7,2) / diffusionLayerDiffusionCoefficient)))
                    #f.write('  Crossing time: %s s\n' % totalCrossingTime)

if __name__ == '__main__':
    main(sys.argv)
