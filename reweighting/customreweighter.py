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
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class CustomReweighter(AbstractReweighter):

    def __init__(self, evaluator):
        ### initializer
        super().__init__()
        self.evaluator = evaluator
        self.evaluator_args = inspect.signature(evaluator).parameters

    def get_weights(self, events, systematic, **kwargs):
        ### internal helper function
        return self.evaluator(events, **kwargs)[systematic]

    def weights(self, events, **kwargs):
        ### get nominal weights
        # (overriding abstract method)
        return self.get_weights(events, 'nominal', **kwargs)

    def weightsup(self, events, untype=None, **kwargs):
        ### get up-varied weights
        # (overriding abstract method)
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'up', **kwargs)

    def weightsdown(self, events, unctype=None, **kwargs):
        ### get down-varied weights
        # (overriding abstract method)
        # note: unctype argument is needed for consistent syntax, but not used
        return self.get_weights(events, 'down', **kwargs)

    def weightsvar(self, events, variation, **kwargs):
        ### get up-varied weights
        # (overriding abstract method)
        # note: variation must be either 'up' or 'down'
        self.check_variation(variation)
        return self.get_weights(events, variation, **kwargs)
