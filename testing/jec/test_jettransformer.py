##############################################
# Test jet transformer for JEC uncertainties #
##############################################

# imports
import sys
import os
import argparse
import numpy as np
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.jetselection import jetselection
from samples.sample import year_from_sample_name
from preprocessing.jettransformer import get_jec_files
from preprocessing.jettransformer import JetTransformer

# input arguments:
parser = argparse.ArgumentParser(description='Test jet transformer for JEC uncertainties')
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
print('Number of events: {}'.format(len(jets)))
print('Number of selected jets: {}'.format(ak.sum(jet_mask)))

# make a jet transformer
jecfiles = get_jec_files(year, unctype='single')
transform = JetTransformer(jecfiles['jec'], jecfiles['junc'])

# make corrected jets
corrected_jets = transform.make_corrected_jets(events, jets=jets)

# print available fields in corrected jets
print(corrected_jets.fields)

# print pt of original and corrected jets
print(ak.flatten(jets.pt))
print(ak.flatten(corrected_jets.pt))
print(ak.flatten(corrected_jets.JES_jes.up.pt))
print(ak.flatten(corrected_jets.JES_jes.down.pt))
