#######################
# Combined reweighter #
#######################

import sys
import os
import inspect
import numpy as np
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


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
        if not isinstance(reweighter, AbstractReweighter):
            msg = 'ERROR: reweighter with name {} (type {})'.format(name, type(reweighter))
            msg += ' does not seem to inherit from AbstractReweighter.'
            raise Exception(msg)
        self.reweighters[name] = reweighter

    def remove_reweighter(self, name):
        self.error_if_not_has_reweighter(name)
        _ = self.reweighters.pop(name)

    # getting reweighters and uncertainties

    def get_reweighter(self, name):
        self.error_if_not_has_reweighter(name)
        return self.reweighters[name]

    def get_unctypes(self, name=None):
        ### get all uncertainty types
        # input arguments:
        # - name: name of a reweighter;
        #         if it is specified, a list with its unctypes is returned,
        #         else a dict matching reweighter names to uncertainty types is returned.
        # note: the list of uncertainty types is not simply equal to the list of reweighters,
        # since some reweighters might contain multiple uncertainties
        # (e.g. statistical and systematic uncertainties)
        uncertainties = {}
        if name is not None: return self.get_reweighter(name).get_unctypes()
        for name, reweighter in self.reweighters.items():
            uncertainties[name] = reweighter.get_unctypes()
        return uncertainties

    def get_variations(self, name=None):
        ### get all variations
        # input arguments:
        # - name: name of a reweighter;
        #         if it is specified, a list with its variations is returned,
        #         else a dict matching reweighter names to variations is returned.
        variations = {}
        if name is not None: return self.get_reweighter(name).get_variations()
        for name, reweighter in self.reweighters.items():
            variations[name] = reweighter.get_variations()
        return variations

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

    # getting nominal and varied weights

    def weights(self, events, **kwargs):
        ### get total nominal per-event weights
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return np.array(weights)

    def weightsup(self, events, reweightername, unctype=None, **kwargs):
        ### get total per-event weights with one reweighter varied up
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            if name==reweightername:
              weights = np.multiply(weights, reweighter.weightsup(events, unctype=unctype, **thiskwargs))
            else: weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return np.array(weights)

    def weightsdown(self, events, reweightername, unctype=None, **kwargs):
        ### get total per-event weights with one reweighter varied down
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            if name==reweightername:
                weights = np.multiply(weights, reweighter.weightsdown(events, unctype=unctype, **thiskwargs))
            else: weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return np.array(weights)

    def weightsvar(self, events, reweightername, variation, **kwargs):
        ### get total per-event weights with one reweighter varied
        weights = np.ones(len(events))
        for name, reweighter in self.reweighters.items():
            thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
            if name==reweightername:
              weights = np.multiply(weights, reweighter.weightsvar(events, variation, **thiskwargs))
            else: weights = np.multiply(weights, reweighter.weights(events, **thiskwargs))
        return np.array(weights)

    def singleweights(self, events, reweightername, **kwargs):
        ### get nominal per-event weights for a single reweighter
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weights(events, **thiskwargs)
        return np.array(weights)

    def singleweightsup(self, events, reweightername, unctype=None, **kwargs):
        ### get per-event weights for a single reweighter varied up
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weightsup(events, unctype=unctype, **thiskwargs)
        return np.array(weights)

    def singleweightsdown(self, events, reweightername, unctype=None, **kwargs):
        ### get per-event weights for a single reweighter varied down
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weightsdown(events, unctype=unctype, **thiskwargs)
        return np.array(weights)

    def singleweightsvar(self, events, reweightername, variation, **kwargs):
        ### get per-event weights for a single reweighter varied
        reweighter = self.reweighters[reweightername]
        thiskwargs = self.select_kwargs(reweighter.weights, kwargs)
        weights = reweighter.weightsvar(events, variation, **thiskwargs)
        return np.array(weights)

    # convenience function for getting nominal weights and all variations

    def allweights(self, events, reweighternames=None, wtype='total', verbose=False, **kwargs):
        ### get nominal weights and all variations
        # input arguments:
        # - events and **kwargs are passed through to the weight functions
        # - reweighternames: list of reweighters to evaluate the variations from (default: all)
        # - wtype: choose from 'total' or 'individual'
        #          in case of 'total', the varied weights are total weights,
        #          where all reweighters are evaluated nominally,
        #          except for one reweighter wich is evaluated with a variation;
        #          in case of 'individual', the varied weights are individual weights,
        #          where the other reweighters have not been taken into account.
        # returns a dict mapping weight set names to np arrays with per-event weights.
        # note: the keys in the returned dict are e.g. 'nominal', 'electronid_statup', etc.
        #       in case of wtype=='individual', the nominal individual weights are also added,
        #       (e.g. 'electronid') as they are needed to calculate the total weights.
        res = {}
        variations = self.get_variations()
        # check arguments
        if wtype not in ['total', 'individual']:
            msg = 'ERROR: wtype {} not recognized'.format(wtype)
            raise Exception(msg)
        if verbose: print('Calculating nominal weights'); sys.stdout.flush()
        res['nominal'] = self.weights(events, **kwargs)
        if verbose: print('Calculating weights for sytematic variations'); sys.stdout.flush()
        # loop over reweighters
        for name in self.reweighters.keys():
            if( reweighternames is not None and name not in reweighternames ): continue
            if verbose: print('  - {}'.format(name)); sys.stdout.flush()
            # calculate individual nominal weight of this reweighter
            weights_single_nominal = self.singleweights(events, name, **kwargs)
            # in case of wtype individual, add individual nominal weights to output
            if wtype=='individual': res[name+'_nom'] = weights_single_nominal
            # calculate total weight without this reweighter
            zero_inds = np.nonzero(weights_single_nominal==0)
            weights_single_nominal[zero_inds] = 1.
            weights_withoutsingle = np.divide(res['nominal'], weights_single_nominal)
            # loop over variations
            for variation in variations[name]:
                if verbose: print('    - {}'.format(variation)); sys.stdout.flush()
                key = '_'.join([name, str(variation)])
                # calculate individual varied weight of this reweighter
                weights_single_variation = self.singleweightsvar(events, name, variation, **kwargs)
                # in case of wtype total, calcualte total varied weight
                if wtype=='total':
                    weights_single_variation = np.multiply(weights_withoutsingle, weights_single_variation)
                # add to output
                res[key] = weights_single_variation
        return res
