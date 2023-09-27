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

class MuonIDReweighter(object):

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

    def get_unctypes(self):
        return self.unctypes
        
    def check_unctype(self, unctype):
        ### internal helper function: check validity of provided uncertainty type
        if not unctype in self.unctypes:
            msg = 'ERROR: uncertainty {} not recognized;'.format(unctype)
            msg += ' allowed values are {}'.format(self.unctypes)
            raise Exception(msg)

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

    def weightsup(self, events, muon_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        if muon_mask is None:
            msg = 'ERROR: muon_mask argument is required.'
            raise Exception(msg)
        muons = events.Muon[muon_mask]
        unc = self.get_muon_uncertainty(muons, unctype=unctype)
        weights = self.get_muon_weights(muons, 'nominal') + unc
        weights = ak.prod(weights, axis=1)
        return weights

    def weightsdown(self, events, muon_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       muon_mask is a keyword argument,
        #       but in practice it is a required argument
        if muon_mask is None:
            msg = 'ERROR: muon_mask argument is required.'
            raise Exception(msg)
        muons = events.Muon[muon_mask]
        unc = self.get_muon_uncertainty(muons, unctype=unctype)
        weights = self.get_muon_weights(muons, 'nominal') - unc
        weights = ak.prod(weights, axis=1)
        return weights
