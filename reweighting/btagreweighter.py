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


class BTagReweighter(object):

    def __init__(self, sffile, normalize=False):
        ### initializer
        self.evaluator = CorrectionSet.from_file(sffile)
        self.jsonmap = 'deepJet_shape'
        self.unctypes_udsgbjets = ['lf', 'hf', 'lfstats1', 'lfstats2', 'hfstats1', 'hfstats2']
        self.unctypes_cjets = ['cferr1', 'cferr2']
        self.unctypes = self.unctypes_udsgbjets + self.unctypes_cjets
        self.normalize = normalize
        self.normalization = None

    def get_jet_weights(self, jets, systematic):
        ### internal helper function: get per-jet weights for specified systematic
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

    def get_unctypes(self):
        return self.unctypes

    def get_jet_weights_rss(self, jets, systematics):
        ### internal helper function: get per-jet weights for root-sum-square of systematics
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
        weights = self.get_jet_weights(jets, systematic)
        weights = ak.prod(weights, axis=1)
        return weights

    def weights(self, events, jet_mask=None):
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

    def weightsup(self, events, jet_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        if jet_mask is None:
            msg = 'ERROR: jet_mask argument is required.'
            raise Exception(msg)
        jets = events.Jet[jet_mask]
        if unctype is not None:
            weights = self.get_event_weights(jets, 'up_'+unctype)
            if self.normalize: weights = self.normalizeweights(weights, ak.num(jets), 'up_'+unctype)
        else:
            systematics = ['up_'+unc for unc in self.unctypes]
            weights = (self.get_jet_weights(jets, 'central') 
                       + self.get_jet_weights_rss(jets, systematics))
            weights = ak.prod(weights, axis=1)
        return weights

    def weightsdown(self, events, jet_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       jet_mask is a keyword argument,
        #       but in practice it is a required argument
        if jet_mask is None:
            msg = 'ERROR: jet_mask argument is required.'
            raise Exception(msg)
        jets = events.Jet[jet_mask]
        if unctype is not None:
            weights = self.get_event_weights(jets, 'down_'+unctype)
            if self.normalize: weights = self.normalizeweights(weights, ak.num(jets), 'down_'+unctype)
        else:
            systematics = ['down_'+unc for unc in self.unctypes]
            weights = (self.get_jet_weights(jets, 'central') 
                       - self.get_jet_weights_rss(jets, systematics))
            weights = ak.prod(weights, axis=1)
        return weights

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
