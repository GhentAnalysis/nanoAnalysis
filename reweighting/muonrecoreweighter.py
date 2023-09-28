############################
# Muon RECO reweighter #
############################
# More information:
# - https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016#RECO_efficiency
# - https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/blob/master/Run2/UL/2016_preVFP/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt_schemaV1.json

# Does not work yet because files provided by Muon POG are in correctionlib schema v1,
# while the current approach seems to work for v2 only.


import sys
import os
import array
import numpy as np
import awkward as ak
from correctionlib._core import CorrectionSet
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class MuonRecoReweighter(AbstractReweighter):

    def __init__(self, sffile):
        ### initializer
        # input arguments:
        # - sffile: path to json file holding the scale factors
        # - year: data-taking year (string format)
        super().__init__()
        self.evaluator = CorrectionSet.from_file(sffile)
        self.jsonmap = 'NUM_TrackerMuons_DEN_genTracks'
        self.unctypes = ['syst', 'stat']
        self.variations = []
        for unctype in self.unctypes:
            self.variations.append(unctype+'_up')
            self.variations.append(unctype+'_down')

    def get_muon_weights(self, muons, valuetype):
        ### internal helper function: retrieve per-muon weights
        # note: correctionlib does not seem to handle jagged arrays,
        #       so need to flatten and unflatten as intermediary steps
        # note: also flattened arrays do not seem to work, not clear why...
        #       use a for-loop for now...
        muons_flat, muons_shape = ak.flatten(muons), ak.num(muons)
        muons_pt = np.array(muons_flat.pt).astype(float)
        muons_eta = np.array(muons_flat.eta).astype(float)
        weights = np.zeros(len(muons_pt))
        for i in range(len(muons_pt)):
            weights[i] = self.evaluator[self.jsonmap].evaluate(
                float(abs(muons_eta[i])), float(muons_pt[i]), valuetype)
        weights = ak.unflatten(weights, muons_shape)
        return weights

    def get_muon_uncertainty_rss(self, muons):
        ### internal helper function: take per-muon root-sum-square of uncertainties
        unc = ak.power(self.get_muon_weights(muons, self.unctypes[0]), 2)
        for unctype in self.unctypes[1:]:
            unc += ak.power(self.get_muon_weights(muons, unctype), 2)
        return ak.sqrt(unc, 0.5)

    def get_muon_uncertainty(self, muons, unctype=None):
        ### internal helper function: get per-muon uncertainty
        if unctype is None: return self.get_muon_uncertainty_rss(muons)
        self.check_unctype(unctype)
        return self.get_muon_weights(muons, unctype)

    def weights(self, events, muon_mask=None):
        ### get nominal per-event weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        if muon_mask is None:
            msg = 'ERROR: muon_mask argument is required.'
            raise Exception(msg)
        muons = events.Muon[muon_mask]
        weights = self.get_muon_weights(muons, 'value')
        weights = ak.prod(weights, axis=1)
        return weights

    def weightsupdown(self, events, upordown, muon_mask=None, unctype=None):
        ### internal helper function to group common code for weightsup and weightsdown
        if muon_mask is None:
            msg = 'ERROR: muon_mask argument is required.'
            raise Exception(msg)
        muons = events.Muon[muon_mask]
        unc = self.get_muon_uncertainty(muons, unctype=unctype)
        sumfactor = 1 if upordown=='up' else -1
        weights = self.get_muon_weights(muons, 'value') + sumfactor*unc
        weights = ak.prod(weights, axis=1)
        return weights

    def weightsup(self, events, muon_mask=None, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'up', muon_mask=muon_mask, unctype=unctype)

    def weightsdown(self, events, muon_mask=None, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'down', muon_mask=muon_mask, unctype=unctype)

    def weightsvar(self, events, variation, muon_mask=None):
        ### get varied per-event weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        # note: variation must be chosen from valid variations,
        #       e.g. 'syst_down', 'stat_up', etc.
        self.check_variation(variation)
        (unctype, upordown) = variation.split('_')
        return self.weightsupdown(events, variation, muon_mask=muon_mask, unctype=unctype)
