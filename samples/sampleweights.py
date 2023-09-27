##############################################
# Tools for reading sample generator weights #
##############################################
# Available info in the "Runs" tree in NanoAOD samples, see here:
# https://cms-nanoaod-integration.web.cern.ch/autoDoc/

# import python modules
import os
import sys
import numpy as np
import uproot
from pathlib import Path

# Note: for now only the nominal weight is stored;
#       to be extended with varied weights later.

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
        if not os.path.exists(samplepath):
            msg = 'ERROR in SampleWeights.__init__:'
            msg += ' file {} does not seem to exist.'.format(samplepath)
            raise Exception(msg)
        with uproot.open(samplepath) as inputfile:
            keys = [key.split(';')[0] for key in inputfile.keys()]
            if 'Runs' not in keys:
                msg = 'ERROR in SampleWeights.__init__:'
                msg += ' "Runs" tree not found in file {};'.format(samplepath)
                msg += ' found only {}.'.format(inputfile.keys())
                raise Exception(msg)
            runs = inputfile['Runs'].arrays(library="np")
        self.genEventCount = np.sum(runs['genEventCount'])
        self.genEventSumw = np.sum(runs['genEventSumw'])
        self.genEventSumw2 = np.sum(runs['genEventSumw2'])
        # (note: need to take sum since multiple entries might be present
        #  if the input file was merged from several input files)
        self.lumiweight = None
        if( self.xsec is not None and self.lumi is not None ):
            self.lumiweight = self.xsec * self.lumi / self.genEventSumw

    def __str__(self):
        res = '--- SampleWeights ---\n'
        for key,val in self.__dict__.items():
            res += '  - {}: {}\n'.format(key, val)
        res = res.strip('\n')
        return res

    def lumiweight(self, xsec=None, lumi=None):
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
