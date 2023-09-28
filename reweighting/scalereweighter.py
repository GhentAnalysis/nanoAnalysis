######################################################
# Renormalization and factorization scale reweighter #
######################################################
# More info:
# - https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/


import sys
import os
import numpy as np
import awkward as ak
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class ScaleReweighter(AbstractReweighter):

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
        self.unctypes = ['muR', 'muF', 'combined']
        self.variations = []
        for unctype in self.unctypes:
            self.variations.append(unctype+'_up')
            self.variations.append(unctype+'_down')
    
    def get_weights(self, events, index):
        if self.rtype=='norm':
            return np.ones(len(events))*self.sampleweights.scaleweights()[index]
        weights = np.array(events.LHEScaleWeight)[:,index]
        if self.rtype=='acceptance':
            weights = weights / self.sampleweights.scaleweights()[index]
        return weights

    def weights(self, events):
        ### get nominal per-event weights
        # (overriding abstract method)
        # note: no nominal reweighting, so return ones.
        return np.ones(len(events))

    def weightsup(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype is optional for syntax, but in practice it is required
        self.check_unctype(unctype)
        if unctype=='muR': return self.get_weights(events,7)
        elif unctype=='muF': return self.get_weights(events,5)
        elif unctype=='combined': return self.get_weights(events,8)

    def weightsdown(self, events, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype is optional for syntax, but in practice it is required
        self.check_unctype(unctype)
        if unctype=='muR': return self.get_weights(events,1)
        elif unctype=='muF': return self.get_weights(events,3)
        elif unctype=='combined': return self.get_weights(events,0)

    def weightsvar(self, events, variation):
        ### get varied per-event weights
        # (overriding abstract method)
        self.check_variation(variation)
        unctype, upordown = variation.split('_')
        if upordown=='up': return self.weightsup(events, unctype=unctype)
        if upordown=='down': return self.weightsdown(events, unctype=unctype)
