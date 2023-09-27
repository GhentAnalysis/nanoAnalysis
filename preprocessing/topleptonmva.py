#####################################################
# Add TOP lepton MVA scores for electrons and muons #
#####################################################
# Based on this implementation:
# https://github.com/HephyAnalysisSW/Analysis/blob/UL/Tools/python/mvaTOPreader.py
# Notes:
# - This module depends on the presence of the variables jetPtRatio and jetBTagDeepFlavor,
#   which are not stored by default in the nanoAOD files.
#   Hence, the module leptonvariables must be run first before this one.

# imports
import os
import numpy as np
import xgboost as xgb
import awkward as ak

class TopLeptonMvaReader(object):

    def __init__(self, year, version, verbose=False):
        self.year = year
        self.version = version

        # check arguments
        if year not in ['2016PreVFP','2016PostVFP','2017','2018']:
            msg = 'ERROR in TopLeptonMvaReader:'
            msg += ' year {} not recognized.'.format(year)
            raise Exception(msg)
        if version not in ['ULv1', 'ULv2']:
            msg = 'ERROR in TopLeptonMvaReader:'
            msg += ' year {} not recognized.'.format(year)
            raise Exception(msg)

        # define working points
        self.wps = {'ULv1': [0.20, 0.41, 0.64, 0.81], 
                    'ULv2': [0.59, 0.81, 0.90, 0.94] }

        # set directory and file names
        weightdir = os.path.join(os.path.dirname(__file__), '../data/leptonmva/weights')
        weightfile = 'TOP'
        if self.version == 'ULv2': weightfile += 'v2'
        diryear = year.replace('20','')
        if year=='2016PreVFP': diryear = '16APV'
        if year=='2016PostVFP': diryear = '16'
        weightfile += 'UL' + diryear + '_XGB.weights.bin'
        elweightfile = os.path.join(weightdir, 'el_' + weightfile)
        muweightfile = os.path.join(weightdir, 'mu_' + weightfile)

        # do printouts if requested
        if verbose:
            print('INFO: loading TOP lepton MVA weights with following properties:')
            print('  - year: {}'.format(year))
            print('  - version: {}'.format(version))
            print('  - electron weights file: {}'.format(elweightfile))
            print('  - muon weights file: {}'.format(muweightfile))
        
        # check if weight files exist
        for f in [elweightfile, muweightfile]:
            if not os.path.exists(f):
                msg = 'ERROR in TopLeptonMvaReader:'
                msg += ' file {} does not exist.'.format(f)
                raise Exception(msg)

        # load weights
        self.electronmva = xgb.Booster() 
        self.electronmva.load_model(elweightfile)
        self.muonmva = xgb.Booster()
        self.muonmva.load_model(muweightfile)
        
    def set_scores(self, events, name='mvaTOP'):
        ### add the lepton MVA scores
        # input arguments:
        # - events: an object of type NanoEventsArray
        # - name: name of the field to add to events.Electron and events.Muon

        # collect scores for electrons
        electron_features = ([
          events.Electron.pt,
          events.Electron.eta,
          events.Electron.jetNDauCharged,
          events.Electron.miniPFRelIso_chg,
          events.Electron.miniPFRelIso_all - events.Electron.miniPFRelIso_chg,
          events.Electron.jetPtRelv2,
          events.Electron.jetPtRatio,
          events.Electron.pfRelIso03_all,
          events.Electron.jetBTagDeepFlavor,
          events.Electron.sip3d,
          np.log(np.abs(events.Electron.dxy)),
          np.log(np.abs(events.Electron.dz)),
          events.Electron.mvaFall17V2noIso
        ])
        if( self.version=='ULv2' ):
            electron_features.append(events.Electron.lostHits)
        # set None values (e.g. when there is no matched jet) to zero
        electron_features = ak.where(ak.is_none(electron_features, axis=2), 0., electron_features)
        # flatten into a numpy array of shape (ninstances, nfeatures)
        counts = ak.num(electron_features[0])
        electron_features = np.transpose(np.array(ak.flatten(electron_features, axis=2)))
        # make c-contiguous
        electron_features = np.ascontiguousarray(electron_features)
        # call xgboost predictor
        electron_scores = self.electronmva.inplace_predict(electron_features)
        # unflatten into an awkward array
        electron_scores = ak.unflatten(electron_scores, counts)
        # set the scores as an additional field for electrons
        events.Electron = ak.with_field(events.Electron, electron_scores, where=name)
        events['Electron'] = events.Electron

        # calculate score for muons
        muon_features = ([
          events.Muon.pt,
          events.Muon.eta,
          events.Muon.jetNDauCharged,
          events.Muon.miniPFRelIso_chg,
          events.Muon.miniPFRelIso_all - events.Muon.miniPFRelIso_chg,
          events.Muon.jetPtRelv2,
          events.Muon.jetPtRatio,
          events.Muon.pfRelIso03_all,
          events.Muon.jetBTagDeepFlavor,
          events.Muon.sip3d,
          np.log(np.abs(events.Muon.dxy)),
          np.log(np.abs(events.Muon.dz)),
          events.Muon.segmentComp
        ])
        # set None values (e.g. when there is no matched jet) to zero
        muon_features = ak.where(ak.is_none(muon_features, axis=2), 0., muon_features)
        # flatten into a numpy array of shape (ninstances, nfeatures)
        counts = ak.num(muon_features[0])
        muon_features = np.transpose(np.array(ak.flatten(muon_features, axis=2)))
        # make c-contiguous
        muon_features = np.ascontiguousarray(muon_features)
        # call xgboost predictor
        muon_scores = self.muonmva.inplace_predict(muon_features)
        # unflatten into an awkward array
        muon_scores = ak.unflatten(muon_scores, counts)
        # set the scores as an additional field for muons
        events.Muon = ak.with_field(events.Muon, muon_scores, where=name)
        events['Muon'] = events.Muon
