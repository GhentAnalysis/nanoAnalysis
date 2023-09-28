########################
# B-tagging reweighter #
########################
# More info:
# - https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
# This page does however only mention that correctionlib+json is recommended,
# but does not give any hint on how to load or use them...
# I could not find any good tutorial, so most of it comes from guessing
# and looking at the json file contents.

# Note: this reweighter is primarily intended to use the 'reshaping' method.
# See method 1d here:
# https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1d_Event_reweighting_using_discr 
# or more info here:
# https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
# It is not yet clear to what extent it will be applicable for other methods.

# Allowed values for systematics seem to be the following:
# - "central"
# - "up_<sysname>" and "down_<sysname>"
#   for sysname in cferr1, cferr2, lf, hf, lfstats1, lfstats2, hfstats1, hfstats2
# - "up_jes" and "down_jes"
# - "up_jes<JEC>" and "down_jes<JEC>"
#   where JEC can be the name of a jet energy correction uncertainty source,
#   see https://github.com/LukaLambrecht/ewkino/blob/ttW/weights/src/ReweighterBTagShape.cc
#   and the links above.

# Todo: speedups


import os
import sys
import numpy as np
import awkward as ak
from correctionlib._core import CorrectionSet
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class BTagReweighter(AbstractReweighter):

    def __init__(self, sffile, normalize=False):
        ### initializer
        super().__init__()
        self.evaluator = CorrectionSet.from_file(sffile)
        self.jsonmap = 'deepJet_shape'
        self.unctypes_udsgbjets = ['lf', 'hf', 'lfstats1', 'lfstats2', 'hfstats1', 'hfstats2']
        self.unctypes_cjets = ['cferr1', 'cferr2']
        self.unctypes = self.unctypes_udsgbjets + self.unctypes_cjets
        self.variations = []
        for unctype in self.unctypes:
            self.variations.append(unctype+'_up')
            self.variations.append(unctype+'_down')
        self.normalize = normalize
        self.normalization = None

    def var_to_unc(self, variation):
        ### internal helper function for change between naming conventions
        self.check_variation(variation)
        if variation.endswith('_up'): return (variation[:-3], 'up')
        if variation.endswith('_down'): return (variation[:-5], 'down')

    def var_to_sys(self, variation):
        ### internal helper function for change between naming conventions
        self.check_variation(variation)
        if variation.endswith('_up'): return 'up_'+variation[:-3]
        if variation.endswith('_down'): return 'down_'+variation[:-5]

    def get_jet_weights(self, jets, systematic):
        ### internal helper function: get per-jet weights for specified systematic
        # input arguments:
        # - jets: a jet collection, e.g. events.Jet
        # - systematic: name of a valid systematic in the BTV convention,
        #   i.e. 'central' for nominal, 'up_[unctype]' or 'down_[unctype]'.
        # note: correctionlib does not seem to handle jagged arrays,
        #       so need to flatten and unflatten as intermediary steps
        # note: also flattened arrays do not seem to work, not clear why...
        #       use a for-loop for now...
        # note: correctionlib does not handle separate systematics for udsg, c and b-jets,
        #       so need to explicitly mask the jets based on flavor before calling evaluate.
        jets_flat, jets_shape = ak.flatten(jets), ak.num(jets)
        jets_pt = np.array(jets_flat.pt).astype(float)
        jets_abseta = np.array(abs(jets_flat.eta)).astype(float)
        jets_flavor = np.array(jets_flat.hadronFlavour).astype(int)
        jets_discr = np.array(jets_flat.btagDeepFlavB).astype(float)
        weights = np.ones(len(jets_flat))
        unctype = systematic.replace('up_','').replace('down_','')
        # case where systematic can be applied only to jets of specific flavor
        # (and need to use central for other jets!)
        if( unctype in self.unctypes_cjets or unctype in self.unctypes_udsgbjets ):
            if unctype in self.unctypes_cjets:
                jetids = np.nonzero(jets_flavor==4)[0]
            elif unctype in self.unctypes_udsgbjets:
                jetids = np.nonzero(((jets_flavor==0) | (jets_flavor==5)))[0]
            for idx in range(len(jets_flat)):
                thissystematic = systematic
                if idx not in jetids: thissystematic = 'central'
                weights[idx] = self.evaluator[self.jsonmap].evaluate(
                  thissystematic, jets_flavor[idx], jets_abseta[idx], jets_pt[idx], jets_discr[idx])
        # 'regular' case where systematic can be applied to all jets
        else:
            for idx in range(len(jets_flat)):
                weights[idx] = self.evaluator[self.jsonmap].evaluate(
                  systematic, jets_flavor[idx], jets_abseta[idx], jets_pt[idx], jets_discr[idx])
        weights = ak.unflatten(weights, jets_shape)
        return weights

    def get_jet_weights_rss(self, jets, systematics):
        ### internal helper function: get per-jet weights for root-sum-square of systematics
        # input arguments:
        # - jets: a jet collection, e.g. events.Jet
        # - systematic: list of names of valid systematics in the BTV convention,
        #   i.e. 'central' for nominal, 'up_[unctype]' or 'down_[unctype]'.
        jets_shape = ak.num(jets)
        nomweights = self.get_jet_weights(jets, 'central')
        nomweights_flat = ak.flatten(nomweights).to_numpy()
        rss = np.zeros(len(nomweights_flat))
        for sys in systematics:
            sysweights = self.get_jet_weights(jets, sys)
            rss += np.power( ak.flatten(sysweights).to_numpy() - nomweights_flat, 2 )
        rss = np.sqrt(rss)
        rss = ak.unflatten( rss, jets_shape )
        return rss

    def get_event_weights(self, jets, systematic):
        ### internal helper function: get per-event weights for specified systematic
        # (by taking the product over the per-jet weights)
        # input arguments:
        # - jets: a jet collection, e.g. events.Jet
        # - systematic: name of a valid systematic in the BTV convention,
        #   i.e. 'central' for nominal, 'up_[unctype]' or 'down_[unctype]'.
        weights = self.get_jet_weights(jets, systematic)
        weights = ak.prod(weights, axis=1)
        return weights

    def weights(self, events, jet_mask=None):
        ### get nominal per-event weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        if jet_mask is None:
            msg = 'ERROR: jet_mask argument is required.'
            raise Exception(msg)
        jets = events.Jet[jet_mask]
        weights = self.get_event_weights(jets, 'central')
        if self.normalize: weights = self.normalizeweights(weights, ak.num(jets), 'central')
        return weights

    def weightsupdown(self, events, upordown, jet_mask=None, unctype=None):
        ### internal helper function to group common code for weightsup and weightsdown
        if jet_mask is None:
            msg = 'ERROR: jet_mask argument is required.'
            raise Exception(msg)
        jets = events.Jet[jet_mask]
        if unctype is not None:
            self.check_unctype(unctype)
            systematic = upordown+'_'+unctype
            weights = self.get_event_weights(jets, systematic)
            if self.normalize: weights = self.normalizeweights(weights, ak.num(jets), systematic)
        else:
            systematics = [upordown+'_'+unc for unc in self.unctypes]
            sumfactor = 1 if upordown=='up' else -1
            weights = (self.get_jet_weights(jets, 'central')
                       + sumfactor*self.get_jet_weights_rss(jets, systematics))
            weights = ak.prod(weights, axis=1)
        return weights

    def weightsup(self, events, jet_mask=None, unctype=None):
        # get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype must be chosen from valid unctypes,
        #       e.g. 'hf', 'lfstats1', 'cferr2', etc.
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'up', jet_mask=jet_mask, unctype=unctype)

    def weightsdown(self, events, jet_mask=None, unctype=None):
        # get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype must be chosen from valid unctypes,
        #       e.g. 'hf', 'lfstats1', 'cferr2', etc.
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'down', jet_mask=jet_mask, unctype=unctype)

    def weightsvar(self, events, variation, jet_mask=None):
        # get varied per-event weights
        # (overriding abstract method)
        # note: variation must be chosen from valid variations,
        #       e.g. 'hf_up', 'lfstats1_down', 'cferr2_up', etc.
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        self.check_variation(variation)
        (unctype, upordown) = self.var_to_unc(variation)
        return self.weightsupdown(events, upordown, jet_mask=jet_mask, unctype=unctype)

    def set_normalization(self, jets, unctypes=None):
        if unctypes is None: unctypes = []
        if unctypes=='all': unctypes = self.unctypes
        systematics = ['central']
        systematics += ['up_'+unc for unc in unctypes]
        systematics += ['down_'+unc for unc in unctypes]
        self.normalization = {}
        # make default normalization
        for sys in systematics:
            self.normalization[sys] = {}
            self.normalization[sys][0] = (1,0)
        # make collection of event indices per number of jets
        jets_shape = ak.num(jets).to_numpy()
        njets_set = sorted(list(set(jets_shape)))
        njets_inds = {}
        for njets in njets_set:
            njets_inds[njets] = np.nonzero(jets_shape==njets)
        # loop over systematics
        for sys in systematics:
            # get event weights for this systematic
            weights = self.get_event_weights(jets, sys)
            # determine average weight per number of jets
            for njets in njets_set:
                thisweights = weights[njets_inds[njets]]
                self.normalization[sys][njets] = (np.mean(thisweights),len(thisweights))

    def normalizeweights(self, weights, njets, systematic):
        if self.normalization is None:
            raise Exception('ERROR: requested normalization, but reweighter was not normalized.')
        if systematic not in self.normalization.keys():
            raise Exception('ERROR: normalization not set for systematic {}'.format(systematic))
        normfactors = np.zeros(len(weights))
        njetkeys = np.array(list(self.normalization[systematic].keys()))
        njetkeys_reverse = njetkeys[::-1]
        for i, njet in enumerate(njets):
            njetkey = next(key for key in njetkeys_reverse if key<=njet)
            normfactors[i] = self.normalization[systematic][njetkey][0]
        return np.divide(weights, normfactors)

    def __str__(self):
        res = 'BTagReweighter\n'
        if self.normalization is None: res += '  [no normalization]'
        else:
            res += '  normalization:\n'
            for sys in self.normalization.keys():
                res += '    {}:\n'.format(sys)
                for njets in self.normalization[sys].keys():
                    res += '      {} jets: {} ({})\n'.format(
                      njets, self.normalization[sys][njets][0],
                      self.normalization[sys][njets][1])
        res = res.strip('\n')
        return res
