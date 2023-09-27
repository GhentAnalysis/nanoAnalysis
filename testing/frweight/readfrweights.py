##########################
# Read fake rate weights #
##########################

# imports
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from samples.sample import year_from_sample_name
import eventselection.lepton_selection_tools as lst
from tools.readfakeratetools import readfrmapfromfile, fakerateweight

# input arguments:
parser = argparse.ArgumentParser(description='Read fake rate weights')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-m', '--muonfrmapfile', required=True, type=os.path.abspath)
parser.add_argument('-e', '--electronfrmapfile', required=True, type=os.path.abspath)
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

# read fake rate maps
electronfrmap = readfrmapfromfile(args.electronfrmapfile, year, 'electron', verbose=True)
muonfrmap = readfrmapfromfile(args.muonfrmapfile, year, 'muon', verbose=True)

# calculate charge flip weights
frweights = fakerateweight(events, electronfrmap, muonfrmap, electron_mask=None, muon_mask=None)
print(frweights)
print(frweights.type)
