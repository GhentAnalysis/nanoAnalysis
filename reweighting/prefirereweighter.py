######################
# Prefire reweighter #
######################
# More info:
# - https://github.com/LukaLambrecht/ewkino/blob/ttW/weights/src/ReweighterPrefire.cc
# - https://cms-nanoaod-integration.web.cern.ch/autoDoc


import sys
import os
import numpy as np
import awkward as ak
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class PrefireReweighter(AbstractReweighter):

    def __init__(self, prefiretype='all'):
        ### initializer
        super().__init__()
        self.prefiretype = prefiretype
        allowed_types = ['all', 'ecal', 'muon']
        if not prefiretype in allowed_types:
          msg = 'ERROR: prefire type {} is not recognized;'.format(prefiretype)
          msg += ' allowed values are {}'.format(allowed_values)
          raise Exception(msg)
        if prefiretype=='muon':
            self.unctypes = ['Stat', 'Syst']
            self.variations = []
            for unctype in self.unctypes:
                self.variations.append(unctype+'_up')
                self.variations.append(unctype+'_down')

    def get_weights(self, events, shift, unctype=None):
        ### internal helper function
        fieldname = 'L1PreFiringWeight'
        varname = ''
        if self.prefiretype=='ecal': varname = 'ECAL_'
        elif self.prefiretype=='muon':
            varname = 'Muon_'
            if shift!='Nom':
                if unctype is None:
                    msg = 'ERROR: must specify an unctype for muon prefiring.'
                    raise Exception(msg)
                self.check_unctype(unctype)
                varname += unctype
        varname += shift
        return getattr(getattr(events, fieldname), varname)

    def weights(self, events):
        ### get nominal per-event weights
        # (overriding abstract method)
        return self.get_weights(events, 'Nom')

    def weightsup(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype is ignored except for muon prefiring
        return self.get_weights(events, 'Up', unctype=unctype)

    def weightsdown(self, events, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype is ignored except for muon prefiring
        return self.get_weights(events, 'Dn', unctype=unctype)

    def weightsvar(self, events, variation):
        ### get varied per-event weights
        # (overriding abstract method)
        # note: variation must either 'up' or 'down'
        #       for general and ecal prefiring,
        #       and 'Syst_up', 'Stat_down', etc.
        #       for muon prefiring.
        self.check_variation(variation)
        unctype = None
        upordown = variation
        if '_' in variation: unctype, upordown = variation.split('_')
        if upordown=='up': upordown = 'Up'
        if upordown=='down': upordown = 'Dn'
        return self.get_weights(events, upordown, unctype=unctype)
