############################
# Read charge flip weights #
############################

# imports
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from samples.sample import year_from_sample_name
from tools.readchargefliptools import readcfmapfromfile, chargeflipweight

# input arguments:
parser = argparse.ArgumentParser(description='Read charge flip weights')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-m', '--cfmapfile', required=True, type=os.path.abspath)
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

# read charge flip map
cfmap = readcfmapfromfile(args.cfmapfile, year, 'electron', verbose=True)

# calculate charge flip weights
cfweights = chargeflipweight(events, cfmap, electron_mask=None, docorrectionfactor=True)
