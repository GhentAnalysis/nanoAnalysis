#####################
# Pileup reweighter #
#####################
# More info:
# - see the meta-info in the json files for some hints on how to use them...


import sys
import os
import numpy as np
import awkward as ak
from correctionlib._core import CorrectionSet
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class PileupReweighter(AbstractReweighter):

    def __init__(self, sffile, year):
        ### initializer
        super().__init__()
        self.evaluator = CorrectionSet.from_file(sffile)
        yeardict = {
          '2016PreVFP': '16',
          '2016PostVFP': '16',
          '2017': '17',
          '2018': '18'
        }
        self.jsonmap = 'Collisions{}_UltraLegacy_goldenJSON'.format(yeardict[year])
    
    def get_weights(self, events, systematic):
        ### internal helper function
        # note: correctionlib does not seem to handle jagged arrays,
        #       so need to flatten and unflatten as intermediary steps
        # note: also flattened arrays do not seem to work, not clear why...
        #       use a for-loop for now...
        ntrueint = events.Pileup.nTrueInt
        weights = np.ones(len(ntrueint))
        for idx in range(len(ntrueint)):
            weights[idx] = self.evaluator[self.jsonmap].evaluate(ntrueint[idx], systematic)
        return weights

    def weights(self, events):
        ### get nominal per-event weights
        # (overriding abstract method)
        return self.get_weights(events, 'nominal')

    def weightsup(self, events, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'up')

    def weightsdown(self, events, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'down')

    def weightsvar(self, events, variation):
        ### get varied per-event weights
        # (overriding abstract method)
        # note: variation must be either 'up' or 'down'
        self.check_variation(variation)
        return self.get_weights(events, variation)
