##################
# PDF reweighter #
##################
# More info:
# - https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/


import sys
import os
import numpy as np
import awkward as ak
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class PDFReweighter(AbstractReweighter):

    def __init__(self, sampleweights, rtype=None):
        ### initializer
        # input arguments
        # - rtype: choose from 'norm', 'acceptance', 'total'
        # - sampleweights: an object of type SampleWeights
        #   (see samples.sampleweights)
        super().__init__()
        self.rtype = rtype
        if(rtype not in ['norm', 'acceptance', 'total']):
            msg = 'ERROR: rtype {} not recognized.'.format(rtype)
            raise Exception(msg)
        self.sampleweights = sampleweights
        if(sampleweights is None):
            msg = 'ERROR: sampleweights not provided.'
            raise Exception(msg)
        self.unctypes = None
        self.variations = list(range(0, sampleweights.nLHEPdfSumw))
    
    def get_weights(self, events, index):
        if self.rtype=='norm':
            return np.ones(len(events))*self.sampleweights.pdfweights()[index]
        weights = np.array(events.LHEPdfWeight)[:,index]
        if self.rtype=='acceptance':
            weights = weights / self.sampleweights.pdfweights()[index]
        return weights

    def weights(self, events):
        ### get nominal per-event weights
        # (overriding abstract method)
        # note: no nominal reweighting, so return ones.
        return np.ones(len(events))

    def weightsup(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: only well-defined in case rtype norm;
        #       for other rtypes there are only variations but no up or down.
        # note: unctype is needed for consistent syntax but not used.
        if self.rtype!='norm':
            msg = 'ERROR: weightsup is not well defined for pdf variations'
            raise Exception(msg)
        maxweight = self.sampleweights.pdfweights(returntype='minmax')[1]
        return np.ones(len(events))*maxweight

    def weightsdown(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: only well-defined in case rtype norm;
        #       for other rtypes there are only variations but no up or down.
        # note: unctype is needed for consistent syntax but not used.
        if self.rtype!='norm':
            msg = 'ERROR: weightsdown is not well defined for pdf variations'
            raise Exception(msg)
        minweight = self.sampleweights.pdfweights(returntype='minmax')[0]
        return np.ones(len(events))*minweight

    def weightsvar(self, events, variation):
        ### get varied per-events weights
        # (overriding abstract method)
        # note: variation must be an integer between 0 and the number of pdf weights
        #       (as given by nLHEPdfSumw)
        variation = int(variation)
        self.check_variation(variation)
        return self.get_weights(events, variation)
