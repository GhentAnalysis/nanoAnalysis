#############################################################################
# Testing script for reading lepton MVA weight files and calculating scores #
#############################################################################

# imports
import sys
import os
from pathlib import Path
import argparse
import numpy as np
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# local imports
sys.path.append(str(Path(__file__).parents[2]))
from preprocessing.topleptonmva import TopLeptonMvaReader
import preprocessing.leptonvariables as lepvars
from samples.sample import year_from_sample_name


# input arguments:
parser = argparse.ArgumentParser(description='Test script for lepton MVA')
parser.add_argument('-i', '--inputfile', required=True)
parser.add_argument('-v', '--version', choices=['ULv1', 'ULv2'], default='ULv1')
parser.add_argument('--entry_start', type=int, default=0)
parser.add_argument('--entry_stop', type=int, default=10)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
year = year_from_sample_name(args.inputfile)
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_start=args.entry_start,
    entry_stop=args.entry_stop,
    schemaclass=NanoAODSchema,
    metadata={'year': year}
).events()

# calculate additional needed input variables
lepvars.add_electron_variables(events, variables=['jetPtRatio','jetBTagDeepFlavor'])
lepvars.add_muon_variables(events, variables=['jetPtRatio','jetBTagDeepFlavor'])

# print out all input variables
# (to be used e.g. for synchronization)
print('Input variables for electrons:')
print(events.Electron.pt)
print(events.Electron.eta)
print(events.Electron.jetNDauCharged)
print(events.Electron.miniPFRelIso_chg)
print(events.Electron.miniPFRelIso_all - events.Electron.miniPFRelIso_chg)
print(events.Electron.jetPtRelv2)
print(events.Electron.jetPtRatio)
print(events.Electron.pfRelIso03_all)
print(events.Electron.jetBTagDeepFlavor)
print(events.Electron.sip3d)
print(np.log(np.abs(events.Electron.dxy)))
print(np.log(np.abs(events.Electron.dz)))
print(events.Electron.mvaFall17V2noIso)
print('Input variables for muons:')
print(events.Muon.pt)
print(events.Muon.eta)
print(events.Muon.jetNDauCharged)
print(events.Muon.miniPFRelIso_chg)
print(events.Muon.miniPFRelIso_all - events.Muon.miniPFRelIso_chg)
print(events.Muon.jetPtRelv2)
print(events.Muon.jetPtRatio)
print(events.Muon.pfRelIso03_all)
print(events.Muon.jetBTagDeepFlavor)
print(events.Muon.sip3d)
print(np.log(np.abs(events.Muon.dxy)))
print(np.log(np.abs(events.Muon.dz)))
print(events.Muon.segmentComp)


# calculate the lepton MVA scores
reader = TopLeptonMvaReader(year, args.version, verbose=True)
reader.set_scores(events)

# print the lepton MVA scores
print('Lepton MVA scores for electrons:')
print(events.Electron.mvaTOP)
print('Lepton MVA scores for muons:')
print(events.Muon.mvaTOP)
