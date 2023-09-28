######################################
# Parton shower (ISR/FSR) reweighter #
######################################
# More info:
# - https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/


import sys
import os
import numpy as np
import awkward as ak
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class PSReweighter(AbstractReweighter):

    def __init__(self):
        ### initializer
        super().__init__()
        self.unctypes = ['isr', 'fsr']
        self.variations = []
        for unctype in self.unctypes:
            self.variations.append(unctype+'_up')
            self.variations.append(unctype+'_down')
         
    def weights(self, events):
        ### get nominal per-event weights
        # (overriding abstract method)
        # note: no nominal reweighting, so return ones
        return np.ones(len(events))

    def weightsup(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype is optional for syntax, but in practice it is required
        self.check_unctype(unctype)
        if unctype=='isr': return np.array(events.PSWeight)[:,0]
        elif unctype=='fsr': return np.array(events.PSWeight)[:,1]

    def weightsdown(self, events, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype is optional for syntax, but in practice it is required
        self.check_unctype(unctype)
        if unctype=='isr': return np.array(events.PSWeight)[:,2]
        elif unctype=='fsr': return np.array(events.PSWeight)[:,3]

    def weightsvar(self, events, variation):
        ### get varied per-event weights
        # (overriding abstract method)
        # note: variation could be 'isr_up', 'fsr_down', etc.
        self.check_variation(variation)
        unctype, upordown = variation.split('_')
        if upordown=='up': return self.weightsup(events, unctype=unctype)
        if upordown=='down': return self.weightsdown(events, unctype=unctype)
