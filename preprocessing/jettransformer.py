##################################################
# Apply jet energy corrections and uncertainties #
##################################################
# Implementation based on this example:
# https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Applying-energy-scale-transformations-with-jetmet_tools

# Note: does not work yet.
# Instead, the best approach is probably to use the official
# tools from JetMET in the skimming step to add additional branches.


import os
import sys
import numpy as np
import awkward as ak
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import FactorizedJetCorrector
from coffea.jetmet_tools import JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack
from coffea.jetmet_tools import CorrectedJetsFactory


def get_jec_files(year, jettype='chs', unctype='single'):
    ### get the correct jec files based on year and other settings
    # note: files and recommendations from here:
    #       https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_MC
    # note: JECs in NanoAOD are already applied,
    #       see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets,
    #       so no nominal corrections but only uncertainties are needed.
    if year=='2016PreVFP': prefix = 'Summer19UL16APV_V7_MC'
    elif year=='2016PostVFP': prefix = 'Summer19UL16_V7_MC'
    elif year=='2017': prefix = 'Summer19UL17_V5_MC'  
    elif year=='2018': prefix = 'Summer19UL18_V5_MC'
    else:
        raise Exception('ERROR: year {} not recognized'.format(year))
    if jettype not in ['chs', 'puppi']:
        raise Exception('ERROR: jettype {} not recognized'.format(jettype))
    if unctype not in ['single', 'split', 'grouped']:
        raise Exception('ERROR: unctype {} not recognized'.format(unctype))
    # base file name
    juncfile = '{}_Uncertainty_AK4PFchs.txt'.format(prefix)
    # modify for uncertainty type
    if unctype in ['split', 'grouped']:
        juncfile = juncfile.replace('Uncertainty', 'UncertaintySources')
        if unctype=='grouped':
            juncfile = 'RegroupedV2_' + juncfile
    # modify for jet type
    if jettype=='puppi':
        juncfile = juncfile.replace('PFchs', 'PFPuppi')
    # make absolute path
    juncfile = os.path.join(os.path.dirname(__file__), '../data/jec/{}/{}'.format(prefix,juncfile))
    return {'jec': [], 'junc': [juncfile]}

 
class JetTransformer(object):

    def __init__(self, jecfiles, juncfiles):
        # empty initialization
        self.jec_stack_names = None
        self.junc_stack_names = None
        self.jec_stack = None
        self.evaluator = None
        self.name_map = None
        # make the jec stack
        self.make_jec_stack(jecfiles, juncfiles)
        # make the name map
        name_map = self.jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetA'] = 'area'
        name_map['ptGenJet'] = 'pt_gen'
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'
        self.name_map = name_map

    def make_jec_stack(self, jecfiles, juncfiles):
        # make a suitable name for each provided jec file
        jec_stack_names = []
        for jecfile in jecfiles:
            name = os.path.basename(jecfile)
            name = name.replace('.jec.txt','').replace('.txt','')
            jec_stack_names.append(name)
        junc_stack_names = []
        for juncfile in juncfiles:
            name = os.path.basename(juncfile)
            name = name.replace('.junc.txt','').replace('.txt','')
            junc_stack_names.append(name)
        self.jec_stack_names = jec_stack_names
        self.junc_stack_names = junc_stack_names
        # make a suitable extractor argument for each provided jec file
        jec_extractor_args = ['* * {}'.format(f) for f in jecfiles+juncfiles]
        # make an extractor
        ext = extractor()
        ext.add_weight_sets(jec_extractor_args)
        ext.finalize()
        # make evaluator
        self.evaluator = ext.make_evaluator()
        # make jec stack
        jec_inputs = ({name: self.evaluator[name] 
          for name in jec_stack_names+junc_stack_names})
        self.jec_stack = JECStack(jec_inputs)

    def make_corrected_jets(self, events, jets=None):
        # prepare jets by adding additional variables
        if jets is None: jets = events.Jet
        jets['pt_raw'] = (1 - jets.rawFactor) * jets.pt
        jets['mass_raw'] = (1 - jets.rawFactor) * jets.mass
        jets['pt_gen'] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
        jets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, jets.pt)[0]
        # make corrector and uncertainty calculator
        events_cache = events.caches[0]
        corrector_inputs = ({name: self.evaluator[name] for name in self.jec_stack_names})
        corrector = FactorizedJetCorrector(**corrector_inputs)
        uncertainties_inputs = ({name: self.evaluator[name] for name in self.junc_stack_names})
        uncertainties = JetCorrectionUncertainty(**uncertainties_inputs)
        jet_factory = CorrectedJetsFactory(self.name_map, self.jec_stack)
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
