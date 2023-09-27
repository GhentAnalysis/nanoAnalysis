######################################################
# Tools for event selection using sample information #
######################################################
# note: largely copied/translated from 
# https://github.com/LukaLambrecht/ewkino/blob/ttW/ttWAnalysis/eventselection/src/eventSelections.cc

# imports
import sys
import os
from pathlib import Path
import awkward as ak
import numpy as np
sys.path.append(str(Path(__file__).parents[1]))
from eventselection.lepton_selection_tools import get_leptons

# internal helper functions
# (just for convenience in the functions below,
#  not meant to be called from outside)

def lepton_from_ME_external_conversion(leptons):
    mask = (
      leptons.matchPdgId==22
      # (to select leptons from external photon conversions?)
      & leptons.isPrompt
      # (to select prompt photons)
      #& leptons.provenanceConversion==0
      # (to select ME photons; not yet implemented)
    )
    return mask

def photon_sample_category(samplename):
    if(
      samplename.startswith('DYJetsToLL') # Z
      or samplename.startswith('TTTo')    # TT
      or samplename.startswith('TTJets')  # TT
      or samplename.startswith('WJets')    # W
    ): return 'inclusive'
    elif(
      samplename.startswith('ZGToLLG')    # Z
      or samplename.startswith('ZGTo2LG') # Z
      or samplename.startswith('TTGamma') # TT
      or samplename.startswith('TTGJets')  # TT
      or samplename.startswith('WGToLNuG') # W
    ): return 'photon'
    else: return 'other'


# external functions
# (meant to be called from outside while doing event selections)

def pass_photon_overlap_removal(events, cat=None, samplename=None, **kwargs):
    ### return a mask for overlap removal between inclusive and photon samples.
    # input arguments:
    # - cat: sample category, either 'inclusive', 'photon' or 'other';
    #   can also be None, in which case the category is determined from the sample name.
    # - samplename: name of the sample from which events were extracted.
    if cat is None:
        if samplename is None:
            msg = 'ERROR: cat and samplename cannot both be None.'
            raise Exception(msg)
        cat = photon_sample_category(samplename)
    if cat=='other': return (events.event >= 0)
    (electrons, muons) = get_leptons(events, **kwargs)
    mask = ak.concatenate((lepton_from_ME_external_conversion(electrons),
      lepton_from_ME_external_conversion(muons)), axis=1)
    mask = (ak.sum(mask,axis=1)>0)
    if cat=='photon': return mask
    elif cat=='inclusive': return ~mask
    else:
        msg = 'ERROR: cat {} not recognized.'.format(cat)
        raise Exception(msg)
