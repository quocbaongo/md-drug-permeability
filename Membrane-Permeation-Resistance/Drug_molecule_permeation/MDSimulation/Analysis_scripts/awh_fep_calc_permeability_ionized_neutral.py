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
from awh_fep_calc_permeability import loadData, calcPermeationResistance, calcCrossingTime

R = 8.3144598e-3

PROGRESS_BAR_LENGTH = 20

def makeUncertaintiesUfloat(data):

    nRows, nCols = data.shape

    if nCols <= 2:
        return data

    newData = numpy.empty((nRows, 2), dtype=object)

    newData[:, 0] = data[:, 0]
    newData[:, 1] = unumpy.uarray(data[:,1], data[:,2])

    return newData

def averageIonizedAndNeutralPmfs(dataArrayList, beta):

    nRows = dataArrayList[0].shape[0]
    nPmfs = len(dataArrayList)
    data = numpy.empty((nPmfs, nRows), dtype=object)
    boltzmannWeights = numpy.empty((nPmfs, nRows), dtype=object)
    xVals = dataArrayList[0][:,0]
    for i in range(nPmfs):
        data[i, :] = dataArrayList[i][:,1]
        boltzmannWeights[i, :] = unumpy.exp(-dataArrayList[i][:,1] * beta)

    xAndResults = numpy.empty((nRows, 2), dtype=object)
    xAndResults[:, 0] = xVals
    xAndResults[:, 1] = numpy.average(data, axis=0, weights=boltzmannWeights)
    #print("%s %s %s %s %s" % (dataArrayList[0][0,1],dataArrayList[1][0,1], boltzmannWeights[0,0], boltzmannWeights[1,0], xAndResults[0,1]))

    return xAndResults, boltzmannWeights

def combineIonizedAndNeutralPmfs(dataArrayList, beta):

    nRows = dataArrayList[0].shape[0]
    nPmfs = len(dataArrayList)
    boltzmannWeights = numpy.empty((nPmfs, nRows), dtype=object)
    xVals = dataArrayList[0][:,0]
    for i in range(nPmfs):
        boltzmannWeights[i, :] = unumpy.exp(-dataArrayList[i][:,1] * beta)

    xAndResults = numpy.empty((nRows, 2), dtype=object)
    xAndResults[:, 0] = xVals
    xAndResults[:, 1] = numpy.sum(boltzmannWeights, axis=0)
    xAndResults[:, 1] = -unumpy.log(xAndResults[:, 1]) / beta
    #print("%s %s %s %s %s" % (dataArrayList[0][0,1],dataArrayList[1][0,1], boltzmannWeights[0,0], boltzmannWeights[1,0], xAndResults[0,1]))

    return xAndResults, boltzmannWeights

def averageDiffusionCoefficients(diffusionArrayList, weights):

    nRows = diffusionArrayList[0].shape[0]
    nCoeffs = len(diffusionArrayList)
    data = numpy.empty((nCoeffs, nRows), dtype=object)
    xVals = diffusionArrayList[0][:,0]
    for i in range(nCoeffs):
        data[i, :] = diffusionArrayList[i][:,1]

    xAndResults = numpy.empty((nRows, 2), dtype=object)
    xAndResults[:, 0] = xVals
    xAndResults[:, 1] = numpy.average(data, axis=0, weights=weights)

    return xAndResults

def main(argv):

    descriptionStr = 'Calculates the weighted PMF and local diffusion coefficient profile from one PMF of a neutral and one from an ionized molecule and the corresponding two local diffusion coefficient profiles.'

    parser=argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-p', '--pmf', type=str, nargs='+',
                        help='Two XVG file(s) containing the pmfs. Must be in the same order as the files specified by --diffusion.')
    parser.add_argument('-d', '--diffusion', type=str, nargs='+',
                        help='Two XVG files containing the diffusion profiles. Must be in the same order as the files specified by --pmf.')
    parser.add_argument('-c', '--concentration', type=float, nargs='+', default=[1, 1],
                        help = 'Concentration ratio. Must be in the same order as the files specified by --pmf, if specified at all. Default 1 1')
    parser.add_argument('-f', '--multiplicationfactor', type=float, nargs='+', default=[1, 1],
                        help = 'Multiplication factor. Multiply the probability by this amount by shifting the PMF. Applied before --concentration if both are used, but it is recommended to use only one of them. Default 1 1')
    parser.add_argument('-t', '--temperature', nargs='?', type=float, default=303.15,
                        help='The temperature of the system. Default 303.15 K')
    parser.add_argument('-z', '--auto_zero_up', action='store_true',
                        help='Automatically shift the PMFs up so that each minimum (including the reference solvent) is zero (after adjusting according to the concentrationRatio).')
    #parser.add_argument('-zz', '--auto_zero_full', action='store_true',
                        #help='Automatically shift the PMFs so that each minimum is zero (after adjusting according to the concentrationRatio).')
    parser.add_argument('-l', '--layers', nargs='+', type=float, default=[1],
                        help='Calculate the permeability coefficient (and lag time) using this number of layers of the system. This can either be a single integer or a pair of floats, with the second being the SEM.')
    parser.add_argument('-o', '--outputbase', type=str,
                        help='The base name of the files to write PMF (xvg), diffusion profile (xvg) and permeation resistance output.')

    args = parser.parse_args()

    pmfFiles = args.pmf
    diffusionFiles = args.diffusion
    concentrationRatio = args.concentration
    multiplicationFactor = args.multiplicationfactor
    outputBase = args.outputbase
    temp = args.temperature
    beta = 1/(R*temp)
    autoZeroUp = args.auto_zero_up
    #autoZeroFull = args.auto_zero_full

    if not pmfFiles or not diffusionFiles:
        parser.error('Must specify --pmf and --diffusion.')
    if len(pmfFiles) != 2 or len(diffusionFiles) != 2:
        parser.error('Must provide two PMFs (from neutral and charged/ionized) as well as two diffusion profiles (from neutral and charged/ionized).')

    multiplicationFactor = numpy.array(multiplicationFactor)
    multiplicationFactorPmfShift = -numpy.log(multiplicationFactor) / beta

    concentrationFactor = numpy.array(concentrationRatio)
    concentrationFactor = concentrationFactor / (concentrationFactor[0] + concentrationFactor[1])
    concentrationFactorPmfShift = -numpy.log(concentrationFactor) / beta

    finalPmfShift = multiplicationFactorPmfShift + concentrationFactorPmfShift

    nLayers = args.layers
    if len(nLayers) == 1:
        nLayers = int(nLayers[0])
    elif len(nLayers) == 2:
        nLayers = ufloat(nLayers[0], nLayers[1])
    else:
        parser.error('nLayers can either be a single integer or a pair of floats where the second one is the SEM of the first one.')

    dataArrayList = []
    for pmfFile in pmfFiles:
        dataArrayList.append(loadData(pmfFile, False))

    pmfMin = ufloat(0,0)

    for i, pmf in enumerate(dataArrayList):
        pmf = makeUncertaintiesUfloat(pmf)
        dataArrayList[i] = pmf
        dataArrayList[i][:,1] += finalPmfShift[i]
        print('Shifting PMF %d by %f based on the relative concentration of the two states.' % (i, finalPmfShift[i]))
        if autoZeroUp:
            pmfMin = min(pmfMin, numpy.amin(dataArrayList[i][:,1]))

    if autoZeroUp:
        for i in range(len(dataArrayList)):
            uncertainty = pmfMin.std_dev

            # Probability of the lowest value being below 0.
            probabilityFactor = max(1e-12, stats.norm.cdf(0, pmfMin.nominal_value, uncertainty))
            shiftValue = min(0, pmfMin.nominal_value)

            uncertainty *= probabilityFactor
            shiftValue = ufloat(shiftValue, uncertainty)
            print("Shifting PMF %d (up to 0) by: %s" % (i, -shiftValue))
            dataArrayList[i][:,1] -= shiftValue

    dataArray, pmfWeights = combineIonizedAndNeutralPmfs(dataArrayList, beta)
    #print(pmfWeights)

    diffusionArrayList = []
    for diffusionFile in diffusionFiles:
        diffusion = loadData(diffusionFile, False)
        diffusion = makeUncertaintiesUfloat(diffusion)
        diffusionArrayList.append(diffusion)

    diffusionArray = averageDiffusionCoefficients(diffusionArrayList, pmfWeights)

    pmfXVals = dataArray[:,0]
    diffusionXVals = diffusionArray[:, 0]
    pmf = dataArray[:,1]
    diffusionProfile = diffusionArray[:,1]

    #print dataArray

    nRows = pmfXVals.shape[0]

    permeationResistance, permeationResistanceProfile = calcPermeationResistance(pmfXVals, pmf, diffusionProfile, beta)



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
        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(pmf), unumpy.std_devs(pmf)))
        numpy.savetxt(fname='%s_pmf.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

        for i in range(len(dataArrayList)):
            data = numpy.column_stack((pmfXVals, unumpy.nominal_values(dataArrayList[i][:,1]), unumpy.std_devs(dataArrayList[i][:,1])))
            numpy.savetxt(fname='%s_pmf_contrib%d.xvg' % (outputBase, i+1), X=data, fmt='%12e', header=header, comments='')

        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(diffusionProfile), unumpy.std_devs(diffusionProfile)))
        numpy.savetxt(fname='%s_diffusion.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

        data = numpy.column_stack((pmfXVals, unumpy.nominal_values(permeationResistanceProfile), unumpy.std_devs(permeationResistanceProfile)))
        numpy.savetxt(fname='%s_perm_resistance.xvg' % outputBase, X=data, fmt='%12e', header=header, comments='')

    permeationResistanceUnc = permeationResistance.std_dev
    permeationResistance = ufloat(permeationResistance.nominal_value, permeationResistanceUnc)
    p = 3600/permeationResistance
    crossingTime = calcCrossingTime(pmfXVals, pmf, p, beta)

    print('PermeationResistance: %s s/cm (%s h/cm)' % (permeationResistance, permeationResistance / 3600))
    print('P: %s cm/h' % p)
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

    if outputBase:
        with open(outputBase+'.txt', 'w') as f:
            f.write(' '.join(argv) + '\n\n')
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

if __name__ == '__main__':
    main(sys.argv)
