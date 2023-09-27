############################################################
# Tools for reading fake rate maps and calculating weights #
############################################################
# translated from https://github.com/LukaLambrecht/ewkino/blob/ttW/Tools/src/readFakeRateTools.cc
# coffea evaluator from this example: https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Coffea-lookup_tools

import sys
import os
import awkward as ak
from coffea.lookup_tools import extractor

def fakerateweight(events, electronfrmap, muonfrmap, electron_mask=None, muon_mask=None):
    # get selected electrons and muons
    electrons = events.Electron
    if electron_mask is not None: electrons = events.Electron[electron_mask]
    muons = events.Muon
    if muon_mask is not None: muons = events.Muon[muon_mask]
    # calculate weights per lepton
    electronpt = ak.where(electrons.pt<44.9, electrons.pt, 44.9)
    efrvals = electronfrmap['frvalue'](electronpt, abs(electrons.eta))
    eweights = -efrvals / (1 - efrvals)
    muonpt = ak.where(muons.pt<44.9, muons.pt, 44.9)
    mfrvals = muonfrmap['frvalue'](muonpt, abs(muons.eta))
    mweights = -mfrvals / (1 - mfrvals)
    # calculate weights per event
    lepweights = ak.concatenate((eweights, mweights), axis=1)
    weights = ak.fill_none(-ak.prod(lepweights, axis=1, mask_identity=True), 0.)
    # printouts for testing
    doprint = False
    if doprint:
        print('Electron pt:')
        print(electrons.pt)
        print('Electron abs eta:')
        print(abs(electrons.eta))
        print('FR values per electron:')
        print(efrvals)
        print('FR weights per electron:')
        print(eweights)
        print('Muon pt:')
        print(muons.pt)
        print('Muon abs eta:')
        print(abs(muons.eta))
        print('FR values per muon:')
        print(mfrvals)
        print('FR weights per muon:')
        print(mweights)
        print('FR weights per event:')
        print(weights)
    return weights

def readfrmapfromfile(filepath, year, flavour,
    fmt='coffea_evaluator', verbose=False):
    # depends on naming conventions, change as needed!
    histname = 'fakeRate_{}_{}'.format(flavour, year)
    frmap = None
    if fmt=='TH2':
        raise Exception('Not supported.')
    elif fmt=='coffea_evaluator':
        ext = extractor()
        ext.add_weight_sets(['frvalue {} {}'.format(histname, filepath)])
        ext.finalize()
        frmap = ext.make_evaluator()
        if verbose:
            print('INFO in readfrmapfromfile: constructed evaluator:')
            print(frmap['frvalue'])
    if frmap is None:
        msg = 'ERROR in readfrmapfromfile:'
        msg += ' format {} not recognized.'.format(fmt)
        raise Exception(msg)
    return frmap

def readfrmap(directory,
    year=None, flavour=None,
    fmt='coffea_evaluator', verbose=False):
    # depends on naming conventions, change as needed!
    filename = "fakeRateMap_data_" + flavour + "_" + year + "_mT.root"
    filepath = os.path.join(directory, filename)
    return readfrmapfromfile(filepath, year, flavour, fmt=fmt, verbose=verbose)
