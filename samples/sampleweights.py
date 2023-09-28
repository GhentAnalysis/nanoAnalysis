##############################################
# Tools for reading sample generator weights #
##############################################
# Available info in the "Runs" tree in NanoAOD samples, see here:
# https://cms-nanoaod-integration.web.cern.ch/autoDoc/
# For the order of the PDF and Scale variations, see the same webpage
# (check the variables LHEPdfWeight and LHEScaleWeight in "Events" tree)

# import python modules
import os
import sys
import numpy as np
import uproot
from pathlib import Path


class SampleWeights(object):
    ### object holding generator weights for a sample
    
    def __init__(self, samplepath, xsec=None, lumi=None):
        ### initializer
        # input arguments:
        # - samplepath: path to a NanoAOD sample file
        #   (should contain at least the "Runs" tree)
        # - xsec: cross-section of the sample
        #   (optional, can also be provided later)
        # - lumi: luminosity to normalize to
        #   (optional, can also be provided later)
        #   (units must correspond to those of xsec!)
        self.samplepath = samplepath
        self.xsec = xsec # can be None
        self.lumi = lumi # can be None
        # check if sample exists
        if not os.path.exists(samplepath):
            msg = 'ERROR in SampleWeights.__init__:'
            msg += ' file {} does not seem to exist.'.format(samplepath)
            raise Exception(msg)
        # read "Runs" tree
        with uproot.open(samplepath) as inputfile:
            keys = [key.split(';')[0] for key in inputfile.keys()]
            if 'Runs' not in keys:
                msg = 'ERROR in SampleWeights.__init__:'
                msg += ' "Runs" tree not found in file {};'.format(samplepath)
                msg += ' found only {}.'.format(inputfile.keys())
                raise Exception(msg)
            runs = inputfile['Runs'].arrays(library="np")
        # read nominal properties
        # (note: need to take sum since multiple entries might be present
        #  if the input file was merged from several input files)
        self.genEventCount = np.sum(runs['genEventCount'])
        self.genEventSumw = np.sum(runs['genEventSumw'])
        self.genEventSumw2 = np.sum(runs['genEventSumw2'])
        # set nominal lumiweight
        self.lumiweight = None
        if( self.xsec is not None and self.lumi is not None ):
            self.lumiweight = self.xsec * self.lumi / self.genEventSumw
        # read pdf and scale variations
        self.nLHEPdfSumw, self.LHEPdfSumw = self.readvariations(runs, 'LHEPdfSumw')
        self.nLHEScaleSumw, self.LHEScaleSumw = self.readvariations(runs, 'LHEScaleSumw')

    def readvariations(self, runs, key):
        ### internal helper function to read cross section variations
        # check that the key is in runs
        nkey = 'n'+key
        if( key not in runs.keys() ):
            msg = 'ERROR: could not find requested key {},'.format(key)
            msg += ' found {}'.format(runs.keys())
            raise Exception(msg)
        if( nkey not in runs.keys() ):
            msg = 'ERROR: could not find requested key {},'.format(nkey)
            msg += ' found {}'.format(runs.keys())
            raise Exception(msg)
        # check that the number of variations is consistent
        if len(set(runs[nkey])) != 1:
            msg = 'ERROR: found inconsistent number of weights'
            msg += ' for {}: {}'.format(key, runs[nkey])
            raise Exception(msg)
        nvariations = runs[nkey][0]
        # read variations
        variations = runs[key]
        nominal = runs['genEventSumw']
        variations = np.sum(np.multiply(variations, nominal))/self.genEventSumw
        return nvariations, variations

    def __str__(self):
        res = '--- SampleWeights ---\n'
        for key,val in self.__dict__.items():
            res += '  - {}: {}\n'.format(key, val)
        res = res.strip('\n')
        return res

    def lumiweight(self, xsec=None, lumi=None):
        ### get event weight scaled with luminosity and cross-section
        # first case: no function arguments,
        # use only info stored in initializer
        if( xsec is None and lumi is None ):
            if self.lumiweight is None:
                msg = 'ERROR in SampleWeights.lumiweight:'
                msg += ' no valid xsec and/or lumi provided.'
                raise Exception(msg)
            return self.lumiweight
        # second case: use xsec and/or lumi from function arguments
        if xsec is None: xsec = self.xsec
        if xsec is None: raise Exception('ERROR: xsec is None')
        if lumi is None: lumi = self.lumi
        if lumi is None: raise Exception('ERROR: lumi is None')
        return xsec * lumi / self.genEventSumw

    def pdfweights(self, returntype='all'):
        ### get pdf weights
        # input arguments:
        # - returntype:
        #   - "all": return all pdf weights
        #   - "minmax": return a tuple of (minimum value, maximum value)
        if returntype=='all': return self.LHEPdfSumw[:]
        elif returntype=='minmax':
            minval = np.min(self.LHEPdfSumw)
            maxval = np.max(self.LHEPdfSumw)
            return (minval, maxval)
        else:
            msg = 'ERROR: return type {} not recognized.'.format(returntype)
            raise Exception(msg)

    def scaleweights(self, returntype='all'):
        ### get scale weights
        # input arguments:
        # - returntype:
        #   - "all": return all pdf weights
        #   - "minmax": return a tuple of (minimum value, maximum value)
        #   - "envelope": same as minmax but exclude some nonphysical combinations
        #                 of muR and muF
        if returntype=='all': return self.LHEScaleSumw[:]
        elif returntype=='minmax':
            minval = np.min(self.LHEScaleSumw)
            maxval = np.max(self.LHEScaleSumw)
            return (minval, maxval)
        elif returntype=='envelope':
            combinations = ([
              # independent variations
              {'muR': 1, 'muF': 0.5},
              {'muR': 1, 'muF': 2},
              {'muR': 0.5, 'muF': 1},
              {'muR': 2, 'muF': 1},
              # positively correlated variations
              {'muR': 0.5, 'muF': 0.5},
              {'muR': 2, 'muF': 2},
              # negatively correlated variations: ignore
            ])
            values = np.array([self.scaleweight(**c) for c in combinations])
            minval = np.min(values)
            maxval = np.max(values)
            return (minval, maxval)
        else:
            msg = 'ERROR: return type {} not recognized.'.format(returntype)
            raise Exception(msg)

    def scaleweight(self, muR=1, muF=1):
        ### get a specific scale weight
        # note: order of the weights in NanoAOD can be retrieve here:
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/
        # input arguments:
        # - muR and muF should be 0.5, 1 or 2
        muR = int(muR) # will convert 0.5 to 0
        muF = int(muF) # will convert 0.5 to 0
        index = 3*muR + muF
        return self.LHEScaleSumw[index]
