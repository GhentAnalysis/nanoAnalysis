###########################################
# Testing script for reading trigger info #
###########################################

# imports
import sys
import os
from pathlib import Path
import argparse
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# local imports
sys.path.append(str(Path(__file__).parents[2]))
from preprocessing.triggervariables import add_trigger_variables 
from samples.sample import year_from_sample_name


# input arguments:
parser = argparse.ArgumentParser(description='Test script for trigger variables')
parser.add_argument('-i', '--inputfile', required=True)
parser.add_argument('--entry_start', type=int, default=0)
parser.add_argument('--entry_stop', type=int, default=-1)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
year = year_from_sample_name(args.inputfile)
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_start=args.entry_start if args.entry_start>=0 else None,
    entry_stop=args.entry_stop if args.entry_stop>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year}
).events()

# add additional required variables
triggerdefs = add_trigger_variables(events, year=year, returntriggerdefs=True)

# define inline function for printing
def tostr(arr):
    return ''.join(['{}'.format(int(el)) for el in arr])

# do printouts for checking
for trigger, hlts in triggerdefs.items():
    print(trigger)
    for hlt in hlts:
        print(tostr(events.HLT[hlt]))
    print('-------')
    print(tostr(events.HLT[trigger]))
    print('-------')
