####################################################
# Tools for lepton pair and Z boson reconstruction #
####################################################
# To do: turn this into an auxiliary class to avoid repeating the same calculations

# imports
import sys
import os
import awkward as ak
import numpy as np
sys.path.append(os.path.abspath('../'))
from constants.particlemasses import m_Z


class ZReco(object):

    def __init__(self, events, dozreco=True, halfwindow=10., samesign=False,
                 electron_mask=None, muon_mask=None):
        ### initializer
        # input arguments:
        # - events: NanoEventsArray read with coffea
        # - dozreco: perform Z boson reconstruction
        #   (can set to False for speed, if only basic sign/flavor info of lepton pairs is needed)
        # - halfwindow: half mass window around Z boson mass (in GeV)
        # - samesign: use same-sign lepton pairs (instead of default opposite-sign)
        # - electron_mask: mask to apply to electrons in events
        # - muon_mask: mask to apply to muons in events
        
        # set parameters
        self.halfwindow = halfwindow
        self.samesign = samesign

        # get leptons
        self.electrons = events.Electron
        if electron_mask is not None: self.electrons = events.Electron[electron_mask]
        self.muons = events.Muon
        if muon_mask is not None: self.muons = events.Muon[muon_mask]

        # do the reconstruction
        self.haszreco = False
        if not dozreco: return
        self.reco()

    def reco(self):
        ### internal helper function for reconstruction
        # get Z->ee and Z->mm candidates
        self.ztoee_candidates = self.reco_from_leptons(self.electrons)
        self.ztomm_candidates = self.reco_from_leptons(self.muons)
        self.ztoll_candidates = ak.concatenate( (self.ztoee_candidates, self.ztomm_candidates), axis=1 )
        # get best Z->ll candidate per event
        mass = (self.ztoll_candidates.i0 + self.ztoll_candidates.i1).mass
        massdev = abs(mass - m_Z)
        sorted_indices = ak.argsort(massdev, axis=1)
        self.ztoll_candidates = self.ztoll_candidates[sorted_indices]
        self.ztoll_best_candidate = ak.pad_none(self.ztoll_candidates, 1, clip=True)
        self.ztoll_best_mass = (self.ztoll_best_candidate.i0 + self.ztoll_best_candidate.i1).mass
        self.haszreco = True

    def reco_if_needed(self):
        ### internal helper function for reconstruction
        if not self.haszreco: self.reco()

    def reco_from_leptons(self, leptons):
        ### internal helper function for reconstruction
        pairs = ak.combinations(leptons, 2, fields=['i0', 'i1'])
        sign_mask = (pairs.i0.charge != pairs.i1.charge)
        if self.samesign: sign_mask = ~sign_mask
        mass_mask = ( abs((pairs.i0 + pairs.i1).mass - m_Z) < self.halfwindow )
        candidates = pairs[ (sign_mask & mass_mask) ]
        return candidates

    def has_os_lepton_pair(self):
        charges = ak.concatenate((self.electrons.charge, self.muons.charge), axis=1)
        return (ak.num(charges) >= 2 
                & abs(ak.sum(charges, axis=1)) != ak.num(charges))

    def has_ossf_electron_pair(self):
        return (ak.num(self.electrons.charge) >= 2
                & abs(ak.sum(self.electrons.charge, axis=1)) != ak.num(self.electrons.charge))

    def has_ossf_muon_pair(self):
        return (ak.num(self.muons.charge) >= 2
                & abs(ak.sum(self.muons.charge, axis=1)) != ak.num(self.muons.charge))

    def has_ossf_lepton_pair(self):
        return (self.has_OSSF_electron_pair() | self.has_OSSF_muon_pair())

    def get_ztoee_candidates(self):
        self.reco_if_needed()
        return self.ztoee_candidates

    def get_ztomm_candidates(self):
        self.reco_if_needed()
        return self.ztomm_candidates

    def get_ztoll_candidates(self):
        self.reco_if_needed()
        return self.ztoll_candidates

    def get_ztoll_best_candidate(self):
        self.reco_if_needed()
        return self.ztoll_best_candidate

    def get_ztoll_best_mass(self):
        self.reco_if_needed()
        return self.ztoll_best_mass

    def has_ztoll_candidate(self):
        candidates = self.get_ztoll_candidates()
        return ( ak.num(candidates) > 0 )

    def n_ztoll_candidates(self):
        candidates = self.get_ztoll_candidates()
        return ak.num(candidates)
