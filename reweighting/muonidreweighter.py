######################
# Muon ID reweighter #
######################
# More information:
# - https://github.com/LukaLambrecht/ewkino/blob/d536cd7082f44875ed61a44282d5bac5c12bfe2d/weights/src/ConcreteReweighterFactory.cc#L245
# - https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Coffea-lookup_tools


import sys
import os
import array
import uproot
import numpy as np
import awkward as ak
from coffea.lookup_tools.dense_lookup import dense_lookup
from pathlib import Path
sys.path.append(Path(__file__).parents[1])
from reweighting.abstractreweighter import AbstractReweighter


class MuonIDReweighter(AbstractReweighter):

    def __init__(self, sffile, wp):
        ### initializer
        # input arguments:
        # - sffile: path to ROOT file holding scale factors
        # - wp: working point (open one of the scale factor files
        #   and check naming convention for allowed values)
        # notes:
        # - the bin errors in the nominal histogram are irrelevant and can be ignored.
        # - the histogram "[...]_combined_syst" contains the absolute (systematic) uncertainties 
        #   as bin contents.
        # - the histogram "[]_stat" contains the absolute (statistical) uncertainties as bin errors.
        super().__init__()
        with uproot.open(sffile) as f:
            basename = 'NUM_LeptonMva{}_DEN_TrackerMuons_abseta_pt'.format(wp)
            (nominal, eta_edges, pt_edges) = f[basename].to_numpy()
            syst = f[basename+'_combined_syst'].values()
            stat = f[basename+'_stat'].errors()
        self.lookup = {
          'nominal': dense_lookup(nominal, [eta_edges, pt_edges]),
          'stat': dense_lookup(stat, [eta_edges, pt_edges]),
          'syst': dense_lookup(syst, [eta_edges, pt_edges])
        }
        self.unctypes = ['syst', 'stat']
        self.variations = []
        for unctype in self.unctypes:
            self.variations.append(unctype+'_up')
            self.variations.append(unctype+'_down')

    def get_muon_weights(self, muons, weighttype):
        ### internal helper function: retrieve per-muon weights
        muons_flat, muons_shape = ak.flatten(muons), ak.num(muons)
        muons_pt = np.array(muons_flat.pt).astype(float)
        muons_eta = np.array(muons_flat.eta).astype(float)
        weights = self.lookup[weighttype](muons_eta, muons_pt)
        weights = ak.unflatten(weights, muons_shape)
        return weights

    def get_muon_uncertainty_rss(self, muons):
        ### internal helper function: retrieve per-muon rss of uncertainties
        unc = self.get_muon_weights(muons, self.unctypes[0])
        muons_shape = ak.num(unc)
        unc = np.power(ak.flatten(unc).to_numpy(), 2)
        for unctype in self.unctypes[1:]:
            unc += np.power( ak.flatten(self.get_muon_weights(muons, unctype)).to_numpy(), 2 )
        unc = np.sqrt(unc)
        return ak.unflatten(unc, muons_shape)

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
        weights = self.get_muon_weights(muons, 'nominal')
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
        weights = self.get_muon_weights(muons, 'nominal') + sumfactor*unc
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
