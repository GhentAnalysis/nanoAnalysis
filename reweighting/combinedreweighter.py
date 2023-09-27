#######################
# Combined reweighter #
#######################

import sys
import os
import inspect
import numpy as np


class CombinedReweighter(object):

    def __init__(self):
        self.reweighters = {}

    # checking presence of reweighters

    def has_reweighter(self, name):
        if name in self.reweighters.keys(): return True
        return False

    def error_if_has_reweighter(self, name):
        if self.has_reweighter(name):
            msg = 'ERROR: a reweighter with name {}'.format(name)
            msg += ' is already present in this combined reweighter.'
            raise Exception(msg)

    def error_if_not_has_reweighter(self, name):
        if not self.has_reweighter(name):
            msg = 'ERROR: a reweighter with name {}'.format(name)
            msg += ' is not present in this combined reweighter.'
            raise Exception(msg)
         
    # adding and removing reweighters
    
    def add_reweighter(self, name, reweighter):
        self.error_if_has_reweighter(name)
        self.reweighters[name] = reweighter

    def remove_reweighter(self, name):
        self.error_if_not_has_reweighter(name)
        _ = self.reweighters.pop(name)

    # getting reweighters and uncertainties

    def get_reweighter(self, name):
        self.error_if_not_has_reweighter(name)
        return self.reweighters[name]

    def get_uncertainties(self, name=None):
        # make a list of all uncertainties.
        # this is not simply equal to the list of reweighters,
        # since some reweighters might contain multiple uncertainties
        # (e.g. statistical and systematic uncertainties)
        uncertainties = {}
        if name is not None: return self.get_reweighter(name).get_unctypes()
        for name, reweighter in self.reweighters.items():
            uncertainties[name] = reweighter.get_unctypes()
        return uncertainties

    # checking arguments of individual reweighters in this combined reweighter

    def get_args(self):
        args = {}
        for name, reweighter in self.reweighters.items():
            args[name] = list(inspect.signature(reweighter.weights).parameters.keys())
        return args

    def select_kwargs(self, f, kwargs):
        ### internal helper function
        # find allowed arguments for function f
        f_args = inspect.signature(f).parameters.keys()
        # if f can take any kwargs, return all kwargs
        if 'kwargs' in f_args: return kwargs
        # else select only allowed arguments
        return {key: val for key, val in kwargs.items() if key in f_args}

    def weights(self, events, **kwargs):
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return weights

    def weightsup(self, events, reweightername, unctype=None, **kwargs):
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            if name==reweightername:
              weights = np.multiply(weights, reweighter.weightsup(events, unctype=unctype, **thiskwargs))
            else: weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return weights

    def weightsdown(self, events, reweightername, unctype=None, **kwargs):
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            if name==reweightername:
                weights = np.multiply(weights, reweighter.weightsdown(events, unctype=unctype, **thiskwargs))
            else: weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return weights

    def singleweights(self, events, reweightername, **kwargs):
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weights(events, **thiskwargs)
        return weights

    def singleweightsup(self, events, reweightername, unctype=None, **kwargs):
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weightsup(events, unctype=unctype, **thiskwargs)
        return weights

    def singleweightsdown(self, events, reweightername, unctype=None, **kwargs):
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weightsdown(events, unctype=unctype, **thiskwargs)
        return weights
