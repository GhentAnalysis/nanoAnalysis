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


class PrefireReweighter(object):

    def __init__(self, prefiretype='all'):
        ### initializer
        self.prefiretype = prefiretype
        allowed_types = ['all', 'ecal', 'muon']
        if not prefiretype in allowed_types:
          msg = 'ERROR: prefire type {} is not recognized;'.format(prefiretype)
          msg += ' allowed values are {}'.format(allowed_values)
          raise Exception(msg)
        self.unctypes = None
        if prefiretype=='muon': self.unctypes=['Stat', 'Syst']

    def get_unctypes(self):
        return self.unctypes

    def check_unctype(self, unctype):
        ### internal helper function: check validity of provided uncertainty type
        if not unctype in self.unctypes:
            msg = 'ERROR: uncertainty {} not recognized;'.format(unctype)
            msg += ' allowed values are {}'.format(self.unctypes)
            raise Exception(msg)

    def get_weights(self, events, shift, unctype=None):
        ### internal helper function
        fieldname = 'L1PreFiringWeight'
        varname = ''
        if self.prefiretype=='ecal': varname = 'ECAL_'
        # todo: extend to statistical up and down for muons
        #       (which are not available for other prefiring types)
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
        return self.get_weights(events, 'Nom')

    def weightsup(self, events, unctype=None):
        return self.get_weights(events, 'Up', unctype=unctype)

    def weightsdown(self, events, unctype=None):
        return self.get_weights(events, 'Dn', unctype=unctype)
