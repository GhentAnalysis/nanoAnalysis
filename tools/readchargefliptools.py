##############################################################
# Tools for reading charge flip maps and calculating weights #
##############################################################
# translated from https://github.com/LukaLambrecht/ewkino/blob/ttW/Tools/src/readChargeFlipTools.cc
# coffea evaluator from this example: https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Coffea-lookup_tools

import sys
import os
import awkward as ak
from coffea.lookup_tools import extractor

def chargeflipweight(events, cfmap, electron_mask=None, docorrectionfactor=False):
    # get selected electrons
    electrons = events.Electron
    if electron_mask is not None: electrons = events.Electron[electron_mask]
    # calculate weights per electron
    cfvals = cfmap['cfvalue'](electrons.pt, abs(electrons.eta))
    eweights = cfvals / (1 - cfvals)
    # calculate weights per event
    sumweights = ak.fill_none(ak.sum(eweights, axis=1, mask_identity=True), 0.)
    prodweights = ak.fill_none(ak.prod(eweights, axis=1, mask_identity=True), 0.)
    # note: we do not need to subtract the product if there was only one electron!
    prodweights = ak.where(ak.count(eweights, axis=1)<2, 0, prodweights)
    weights = sumweights - prodweights
    # printouts for testing
    doprint = False
    if doprint:
        for i in range(ak.count(events.event)):
            print('chargeflipweight results for event {}'.format(events.event[i]))
            print('  Electron pt: {}'.format(electrons.pt[i]))
            print('  Electron abs eta: {}'.format(abs(electrons.eta[i])))
            print('  CF values: {}'.format(cfvals[i]))
            print('  CF weights per electron: {}'.format(eweights[i]))
            print('  CF weights per event:')
            print('    sum: {}'.format(sumweights[i]))
            print('    prod: {}'.format(prodweights[i]))
            print('    tot: {}'.format(weights[i]))
    # modify weights by a correction factor
    if docorrectionfactor:
        year = events.metadata['year']
        if( year=="2016PreVFP" ): weights = 0.85 * weights
        elif( year=="2016PostVFP" ): weights = 0.95 * weights
        elif( year=="2017" ): weights = 1.4 * weights
        elif( year=="2018" ): weights = 1.4 * weights
        else:
            msg = 'ERROR in chargeflipweight:'
            msg += ' year {} not recognized (needed for correction factor).'.format(year)
            raise Exception(msg)
    return weights

def readcfmapfromfile(filepath, year, flavour,
    fmt='coffea_evaluator', verbose=False):
    # depends on naming conventions, change as needed!
    histname = 'chargeFlipRate_{}_{}'.format(flavour, year)
    cfmap = None
    if fmt=='TH2':
        raise Exception('Not supported.')
    elif fmt=='coffea_evaluator':
        ext = extractor()
        ext.add_weight_sets(['cfvalue {} {}'.format(histname, filepath)])
        ext.finalize()
        cfmap = ext.make_evaluator()
        if verbose:
            print('INFO in readcfmapfromfile: constructed evaluator:')
            print(cfmap['cfvalue'])
    if cfmap is None:
        msg = 'ERROR in readcfmapfromfile:'
        msg += ' format {} not recognized.'.format(fmt)
        raise Exception(msg)
    return cfmap

def readcfmap(directory,
    year=None, flavour=None, process=None, binning=None,
    fmt='coffea_evaluator', verbose=False):
    # depends on naming conventions, change as needed!
    filename = "chargeFlipMap_MC_" + flavour + "_" + year
    filename += "_process_" + process + "_binning_" + binning + ".root"
    filepath = os.path.join(directory, filename)
    return readcfmapfromfile(filepath, year, flavour, fmt=fmt, verbose=verbose)
