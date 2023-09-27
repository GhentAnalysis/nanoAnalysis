#################################
# Definition of event variables #
#################################


# import python modules
import sys
import os
from pathlib import Path
import awkward as ak
import numpy as np
# import framework modules
sys.path.append(str(Path(__file__).parents[2]))
from tools.readfakeratetools import fakerateweight
from tools.readchargefliptools import chargeflipweight
from eventreconstruction.zreco import ZReco

def calculate_event_variables(events,
    weights=None, nentries_reweight=1,
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    electronfrmap=None, muonfrmap=None,
    electroncfmap=None ):
    ### calculate event variables
    res = {}
    # initializations
    nevents = ak.count(events.event)
    zreco = ZReco(events, halfwindow=10.,
                  electron_mask=electron_fo_mask,
                  muon_mask=muon_fo_mask)
    # event identifiers
    res['event'] = events.event
    res['run'] = events.run
    res['luminosityBlock'] = events.luminosityBlock
    # total yield (fixed arbitrary value)
    res['yield'] = np.ones(nevents)*0.5
    # generator event weights
    if( events.metadata['dtype']=='sim' ):
        res['genWeight'] = events.genWeight * nentries_reweight
        if( weights is None ):
            msg = 'WARNING in calculate_event_variables:'
            msg += ' no valid SampleWeights object found,'
            msg += ' will not write normalized sample weights.'
            print(msg)
        else:
            res['genNormWeight'] = events.genWeight / weights.genEventSumw * nentries_reweight
    # fake rate and chargeflip rate event weights
    res['fakeRateWeight'] = np.ones(nevents)
    if(electronfrmap is not None and muonfrmap is not None):
        electron_fr_mask = (electron_fo_mask & ~electron_tight_mask)
        muon_fr_mask = (muon_fo_mask & ~muon_tight_mask)
        res['fakeRateWeight'] = fakerateweight(events, electronfrmap, muonfrmap,
          electron_mask=electron_fr_mask, muon_mask=muon_fr_mask)
    res['chargeFlipWeight'] = np.ones(nevents)
    if(electroncfmap is not None):
        res['chargeFlipWeight'] = chargeflipweight(events, electroncfmap,
          electron_mask=electron_tight_mask, docorrectionfactor=True)
    # get object collections
    electrons = events.Electron[electron_fo_mask]
    muons = events.Muon[muon_fo_mask]
    jets = events.Jet[jet_mask]
    bjets = events.Jet[bjet_mask]
    # number of jets
    njets = ak.sum(jet_mask, axis=1)
    res['nJets'] = njets
    # jet pt and eta
    jet_pt = np.array(ak.fill_none(ak.pad_none(jets.pt, 2, axis=1, clip=True), 0.))
    jet_eta = np.array(ak.fill_none(ak.pad_none(jets.eta, 2, axis=1, clip=True), 0.))
    res['jetPtLeading'] = jet_pt[:,0]
    res['jetEtaLeading'] = jet_eta[:,0]
    res['jetPtSubLeading'] = jet_pt[:,1]
    res['jetEtaSubLeading'] = jet_eta[:,1]
    # number of b-jets
    nbjets = ak.sum(bjet_mask, axis=1)
    res['nBJets'] = nbjets
    # b-jet pt and eta
    bjet_pt = np.array(ak.fill_none(ak.pad_none(bjets.pt, 2, axis=1, clip=True), 0.))
    bjet_eta = np.array(ak.fill_none(ak.pad_none(bjets.eta, 2, axis=1, clip=True), 0.))
    res['bjetPtLeading'] = bjet_pt[:,0]
    res['bjetEtaLeading'] = bjet_eta[:,0]
    res['bjetPtSubLeading'] = bjet_pt[:,1]
    res['bjetEtaSubLeading'] = bjet_eta[:,1]
    # MET
    res['MET_pt'] = events.MET.pt
    res['MET_phi'] = events.MET.phi
    # number of muons
    res['nMuons'] = ak.sum(muon_fo_mask, axis=1)
    # number of electrons
    res['nElectrons'] = ak.sum(electron_fo_mask, axis=1)
    # scalar pt sums
    res['HT'] = ak.sum(jets.pt, axis=1)
    res['LT'] = (ak.sum(ak.concatenate((electrons.pt, muons.pt), axis=1), axis=1)
                  + events.MET.pt)
    # number of Z boson candidates
    nz = zreco.n_ztoll_candidates()
    # custom categorization variable for jets and b-jets
    njnb = -np.ones(nevents)
    njnb = ak.where(nbjets==0, np.clip(njets,0,4), njnb)
    njnb = ak.where(nbjets==1, 5+np.clip(njets,0,5)-1, njnb)
    njnb = ak.where(nbjets>1, 10+np.clip(njets,0,5)-2, njnb)
    res['nJetsNBJetsCat'] = njnb
    # custom categorization variable for jets, b-jets and Z candidates
    njnz = -np.ones(nevents)
    njnz = np.where(nz==2, 0, njnz)
    njnz = np.where(nz==1, np.clip(njets,0,2)+1, njnz)
    res['nJetsNZCat'] = njnz 
    return res
