##################################################
# Abstract base class for all reweighter classes #
##################################################

from abc import ABC, abstractmethod
import sys
import os


class AbstractReweighter(ABC):

    def __init__(self):
        ### initializer
        # The default behaviour for a reweighter is to have 
        # no specified uncertainty types,
        # corresponding to just two variations (up and down).
        # This can be overridden in the concrete reweighters,
        # e.g. for split statistical and systematic variations.
        self.unctypes = None
        self.variations = ['up', 'down']

    def get_unctypes(self):
        # should not be overridden in child classes, just inherit
        return self.unctypes

    def check_unctype(self, unctype):
        # should not be overridden in child classes, just inherit
        if not unctype in self.unctypes:
            msg = 'ERROR: uncertainty {} not recognized;'.format(unctype)
            msg += ' allowed values are {}'.format(self.unctypes)
            raise Exception(msg)

    def get_variations(self):
        # should not be overridden in child classes, just inherit
        return self.variations

    def check_variation(self, variation):
        # should not be overridden in child classes, just inherit
        if not variation in self.variations:
            msg = 'ERROR: variation {} not recognized;'.format(variation)
            msg += ' allowed values are {}'.format(self.variations)
            raise Exception(msg)

    @abstractmethod
    def weights(self, events, **kwargs):
        pass

    @abstractmethod
    def weightsup(self, events, untype=None, **kwargs):
        pass

    @abstractmethod
    def weightsdown(self, events, unctype=None, **kwargs):
        pass
 
    @abstractmethod
    def weightsvar(self, events, variation, **kwargs):
        pass
