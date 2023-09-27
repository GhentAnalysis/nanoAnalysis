############################################################
# Testing script for reading lepton selection and cleaning #
############################################################

# imports
import sys
import os
from pathlib import Path
import argparse
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# local imports
sys.path.append(str(Path(__file__).parents[2]))
from preprocessing.leptongenvariables import add_electron_gen_variables
from preprocessing.leptongenvariables import add_muon_gen_variables
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons


# input arguments:
parser = argparse.ArgumentParser(description='Test script for lepton ID')
parser.add_argument('-i', '--inputfile', required=True)
parser.add_argument('--entry_start', type=int, default=0)
parser.add_argument('--entry_stop', type=int, default=-1)
parser.add_argument('--electronid', default=None)
parser.add_argument('--muonid', default=None)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_start=args.entry_start if args.entry_start>=0 else None,
    entry_stop=args.entry_stop if args.entry_stop>=0 else None,
    schemaclass=NanoAODSchema,
).events()

# add additional required variables
add_electron_gen_variables(events, variables=['isPrompt'])
add_muon_gen_variables(events, variables=['isPrompt'])

# define masks for chosen lepton IDs
elselmask = electronselection(events.Electron, selectionid=args.electronid)
muselmask = muonselection(events.Muon, selectionid=args.muonid)

# basic printouts
print('Number of events:')
print(ak.count(events.event))
print('Electrons before selection:')
print(ak.num(events.Electron.pt))
print(ak.count(events.Electron.pt))
print('Electrons after selection:')
print(ak.num(events.Electron[elselmask].pt))
print(ak.count(events.Electron[elselmask].pt))
print('Muons before selection:')
print(ak.num(events.Muon.pt))
print(ak.count(events.Muon.pt))
print('Muons after selection:')
print(ak.num(events.Muon[muselmask].pt))
print(ak.count(events.Muon[muselmask].pt))

# define mask for cleaning
elcleanmask = clean_electrons_from_muons(events.Electron, events.Muon[muselmask])

# additional printouts
print('Electrons after cleaning from muons:')
print(ak.num(events.Electron[elselmask & elcleanmask].pt))
print(ak.count(events.Electron[elselmask & elcleanmask].pt))

# define mask for prompt leptons
elpromptmask = events.Electron.isPrompt
mupromptmask = events.Muon.isPrompt

# additional printouts
print('Prompt electrons before selection:')
print(ak.num(events.Electron[elpromptmask].pt))
print(ak.count(events.Electron[elpromptmask].pt))
print('Prompt electrons after selecion and cleaning:')
print(ak.num(events.Electron[elselmask & elcleanmask & elpromptmask].pt))
print(ak.count(events.Electron[elselmask & elcleanmask & elpromptmask].pt))
print('Prompt muons before selection:')
print(ak.num(events.Muon[mupromptmask].pt))
print(ak.count(events.Muon[mupromptmask].pt))
print('Prompt muons after selecion:')
print(ak.num(events.Muon[muselmask & mupromptmask].pt))
print(ak.count(events.Muon[muselmask & mupromptmask].pt))
