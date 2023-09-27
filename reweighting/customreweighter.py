#####################
# Custom reweighter #
#####################
# This reweighter type allows to add reweighting and/or uncertainties
# based on any given variable.
# The typical use case is to add x % of uncertainty for eventvariables
# with more than y jets / b-tagged jets.
 

import sys
import os
import numpy as np
import awkward as ak
import inspect


class CustomReweighter(object):

    def __init__(self, evaluator):
        ### initializer
        self.evaluator = evaluator
        self.evaluator_args = inspect.signature(evaluator).parameters
        self.unctypes = None

    def get_unctypes(self):
        return self.unctypes

    def get_weights(self, events, systematic, **kwargs):
        ### internal helper function
        return self.evaluator(events, **kwargs)[systematic]

    def weights(self, events, **kwargs):
        return self.get_weights(events, 'nominal', **kwargs)

    def weightsup(self, events, untype=None, **kwargs):
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'up', **kwargs)

    def weightsdown(self, events, unctype=None, **kwargs):
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'down', **kwargs)
