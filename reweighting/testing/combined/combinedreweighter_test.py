############################
# Test combined reweighter #
############################

# imports
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[3]))
from objectselection.jetselection import jetselection
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from reweighting.muonidreweighter import MuonIDReweighter
from reweighting.electronidreweighter import ElectronIDReweighter
from reweighting.combinedreweighter import CombinedReweighter


# input arguments:
parser = argparse.ArgumentParser(description='Test custom reweighter')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
print('Loading events from input file...')
year = year_from_sample_name(args.inputfile)
samplename = os.path.basename(args.inputfile)
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year}
).events()
print('Number of events in input file: {}'.format(ak.count(events.event)))

# calculate object masks
print('Performing jet selection...')
jet_mask = jetselection(events.Jet, selectionid='run2ul_default')
jets = events.Jet[jet_mask]
njets = ak.num(jets)
print('Number of events: {}'.format(len(events)))
print('Number of selected jets: {}'.format(ak.sum(jet_mask)))

# make combined reweighter
muonrecofile = '../../data/leptonid/leptonMVAUL_SF_muons_Medium_2018.root'
electronrecofile = '../../data/leptonid/leptonMVAUL_SF_electrons_Medium_2018.root'
reweighter = CombinedReweighter()
reweighter.add_reweighter('muonid', MuonIDReweighter(muonrecofile, 'Medium'))
reweighter.add_reweighter('electronid', ElectronIDReweighter(electronrecofile))

# printouts for testing
print(reweighter.get_uncertainties())
