##########################
# Electron ID reweighter #
##########################
# More information:
# - https://github.com/LukaLambrecht/ewkino/blob/d536cd7082f44875ed61a44282d5bac5c12bfe2d/weights/src/ConcreteReweighterFactory.cc#L288
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


class ElectronIDReweighter(AbstractReweighter):

    def __init__(self, sffile):
        ### initializer
        # input arguments:
        # - sffile: path to ROOT file holding scale factors
        # notes:
        # - the bin errors in the nominal histogram are irrelevant and can be ignored.
        # - the histogram "sys" contains the absolute (systematic) uncertainties as bin contents.
        # - the histogram "stat" contains the absolute (statistical) uncertainties as bin contents.
        super().__init__()
        with uproot.open(sffile) as f:
            (nominal, eta_edges, pt_edges) = f['EGamma_SF2D'].to_numpy()
            syst = f['sys'].values()
            stat = f['stat'].values()
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

    def get_electron_weights(self, electrons, weighttype):
        ### internal helper function: retrieve per-electron weights
        electrons_flat, electrons_shape = ak.flatten(electrons), ak.num(electrons)
        electrons_pt = np.array(electrons_flat.pt).astype(float)
        electrons_eta = np.array(electrons_flat.eta).astype(float)
        weights = self.lookup[weighttype](electrons_eta, electrons_pt)
        weights = ak.unflatten(weights, electrons_shape)
        return weights

    def get_electron_uncertainty_rss(self, electrons):
        ### internal helper function: retrieve per-electron rss of uncertainties
        unc = self.get_electron_weights(electrons, self.unctypes[0])
        electrons_shape = ak.num(unc)
        unc = np.power(ak.flatten(unc).to_numpy(), 2)
        for unctype in self.unctypes[1:]:
            unc += np.power( ak.flatten(self.get_electron_weights(electrons, unctype)).to_numpy(), 2 )
        unc = np.sqrt(unc)
        return ak.unflatten(unc, electrons_shape)

    def get_electron_uncertainty(self, electrons, unctype=None):
        ### internal helper function: get per-electron uncertainty
        if unctype is None: return self.get_electron_uncertainty_rss(electrons)
        self.check_unctype(unctype)
        return self.get_electron_weights(electrons, unctype)

    def weights(self, events, electron_mask=None):
        ### get nominal weights
        # (overriding abstract method)
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        if electron_mask is None:
            msg = 'ERROR: electron_mask argument is required.'
            raise Exception(msg)
        electrons = events.Electron[electron_mask]
        weights = self.get_electron_weights(electrons, 'nominal')
        weights = ak.prod(weights, axis=1)
        return weights

    def weightsupdown(self, events, upordown, electron_mask=None, unctype=None):
        ### internal helper function to group common code for weightsup and weightsdown
        if electron_mask is None:
            msg = 'ERROR: electron_mask argument is required.'
            raise Exception(msg)
        electrons = events.Electron[electron_mask]
        unc = self.get_electron_uncertainty(electrons, unctype=unctype)
        sumfactor = 1 if upordown=='up' else -1
        weights = self.get_electron_weights(electrons, 'nominal') + sumfactor*unc
        weights = ak.prod(weights, axis=1)
        return weights

    def weightsup(self, events, electron_mask=None, unctype=None):
        ### get up-varied per-event weights
        # (overriding abstract method)
        # note: unctype must be either 'syst' or 'stat'
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'up', electron_mask=electron_mask, unctype=unctype)

    def weightsdown(self, events, electron_mask=None, unctype=None):
        ### get down-varied per-event weights
        # (overriding abstract method)
        # note: unctype must be either 'syst' or 'stat'
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        return self.weightsupdown(events, 'down', electron_mask=electron_mask, unctype=unctype)

    def weightsvar(self, events, variation, electron_mask=None):
        # get varied per-event weights
        # (overriding abstract method)
        # note: variation must be chosen from valid variations,
        #       e.g. 'syst_up', 'stat_down', etc.
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        self.check_variation(variation)
        (unctype, upordown) = variation.split('_')
        return self.weightsupdown(events, upordown, electron_mask=electron_mask, unctype=unctype) 
