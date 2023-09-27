#####################################################
# Test JEC uncertainties and related jet operations #
#####################################################

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
from objectselection.jetuncs import get_varied_jets, get_varied_met
from objectselection.jetuncs import get_available_jec_variations
from objectselection.jetuncs import get_available_jer_variations
from samples.sample import year_from_sample_name

# input arguments:
parser = argparse.ArgumentParser(description='Test JEC uncertainties')
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

# print available variations
print(get_available_jec_variations(events))
print(get_available_jer_variations(events))

# get nominal jets
jets_nominal = events.Jet
print('Jet.pt:')
print(jets_nominal.pt)

# get varied jets
variation = 'jesTotal'
jets_up = get_varied_jets(events, variation+'Up')
jets_down = get_varied_jets(events, variation+'Down')
print('Jet pt up variation:')
print(jets_up.pt)
print('Jet pt down variation:')
print(jets_down.pt)

# check if original pt was not modified
print('Check original pt:')
print(jets_nominal.pt)
print(events.Jet.pt)
print('---------\n')

# same checks for mass
print(jets_nominal.mass)
print(jets_up.mass)
print(jets_down.mass)
print('---------\n')

# get nominal met
met_nominal = events.MET
print('MET.pt')
print(met_nominal.pt)

# get varied jets
variation = 'jesTotal'
met_up = get_varied_met(events, variation+'Up')
met_down = get_varied_met(events, variation+'Down')
print('MET up variation:')
print(met_up.pt)
print('MET down variation:')
print(met_down.pt)

# check if original pt was not modified
print('Check original pt:')
print(met_nominal.pt)
print(events.MET.pt)
print('---------\n')

# same checks for phi
print(met_nominal.phi)
print(met_up.phi)
print(met_down.phi)
