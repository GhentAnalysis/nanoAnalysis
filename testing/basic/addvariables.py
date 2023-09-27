##############################################################
# Test script for adding variables to NanoEventsArray object #
##############################################################

# imports
import sys
import os
from pathlib import Path
import argparse
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# local imports
sys.path.append(str(Path(__file__).parents[2]))
import preprocessing.leptonvariables as lepvars


# input arguments:
parser = argparse.ArgumentParser(description='Test script for adding variables')
parser.add_argument('-i', '--inputfile', required=True)
parser.add_argument('-v', '--variables', nargs='+', default=['all'])
parser.add_argument('--entry_start', type=int, default=0)
parser.add_argument('--entry_stop', type=int, default=10)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_start=args.entry_start,
    entry_stop=args.entry_stop,
    schemaclass=NanoAODSchema,
).events()

# add lepton variables
lepvars.add_electron_variables(events, variables=args.variables)
lepvars.add_muon_variables(events, variables=args.variables)

# test jetPtRatio
#print(events.Electron.jetRelIso)
print(events.Electron.jetPtRatio)
#print(events.Muon.jetRelIso)
print(events.Muon.jetPtRatio)

# test jetBTagDeepFlavor
print(events.Electron.jetBTagDeepFlavor)
print(events.Muon.jetBTagDeepFlavor)
