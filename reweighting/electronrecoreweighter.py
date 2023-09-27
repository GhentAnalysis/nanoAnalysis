############################
# Electron RECO reweighter #
############################
# More information:
# - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# - https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON

import sys
import os
import array
import numpy as np
import awkward as ak
from correctionlib._core import CorrectionSet


class ElectronRecoReweighter(object):

    def __init__(self, sffile, year):
        # input arguments:
        # - sffile: path to json file holding the scale factors
        # - year: data-taking year (string format)
        self.evaluator = CorrectionSet.from_file(sffile)
        self.year = year
        if year=='2016PreVFP': self.year = '2016preVFP'
        if year=='2016PostVFP': self.year = '2016postVFP'
        self.jsonmap = 'UL-Electron-ID-SF'
        self.unctypes = None

    def get_unctypes(self):
        return self.unctypes

    def get_electrons_ptbin(self, electrons, ptbin):
        ### internal helper function
        if ptbin == 'RecoBelow20':
           electrons = electrons[(electrons.pt < 20.)]
        elif ptbin == 'RecoAbove20':
           electrons = electrons[(electrons.pt > 20.)]
        else:
           raise Exception('ERROR: ptbin {} not recognized.'.format(ptbin))
        return electrons
        
    def get_weights(self, electrons, valuetype):
        ### internal helper function
        # note: correctionlib does not seem to handle jagged arrays,
        #       so need to flatten and unflatten as intermediary steps
        # note: also flattened arrays do not seem to work, not clear why...
        #       use a for-loop for now...
        weights = np.ones(len(electrons))
        for ptbin in ['RecoBelow20', 'RecoAbove20']:
            electrons_ptbin = self.get_electrons_ptbin(electrons, ptbin)
            electrons_flat, electrons_shape = ak.flatten(electrons_ptbin), ak.num(electrons_ptbin)
            electrons_pt = np.array(electrons_flat.pt).astype(float)
            electrons_eta = np.array(electrons_flat.eta).astype(float)
            thisweights = np.zeros(len(electrons_pt))
            for i in range(len(electrons_pt)):
                thisweights[i] = self.evaluator[self.jsonmap].evaluate(
                  self.year, valuetype, ptbin, float(electrons_eta[i]), float(electrons_pt[i]))
            #thisweights = self.evaluator[self.jsonmap].evaluate(
            #      self.year, valuetype, ptbin, electrons_eta, electrons_pt)
            thisweights = ak.unflatten(thisweights, electrons_shape)
            thisweights = ak.prod(thisweights, axis=1)
            weights = np.multiply(weights, thisweights)
        return weights

    def weights(self, events, electron_mask=None):
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        if electron_mask is None:
            msg = 'ERROR: electron_mask argument is required.'
            raise Exception(msg)
        electrons = events.Electron[electron_mask]
        return self.get_weights(electrons, 'sf')

    def weightsup(self, events, electron_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        # note: unctype argument is needed for consistent syntax, but not used
        if electron_mask is None:
            msg = 'ERROR: electron_mask argument is required.'
            raise Exception(msg)
        electrons = events.Electron[electron_mask]
        return self.get_weights(electrons, 'sfup')

    def weightsdown(self, events, electron_mask=None, unctype=None):
        # note: for consistent syntax across reweighters,
        #       electron_mask is a keyword argument,
        #       but in practice it is a required argument
        # note: unctype argument is needed for consistent syntax, but not used
        if electron_mask is None:
            msg = 'ERROR: electron_mask argument is required.'
            raise Exception(msg)
        electrons = events.Electron[electron_mask]
        return self.get_weights(electrons, 'sfdown')
