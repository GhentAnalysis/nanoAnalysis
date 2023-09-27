#############################
# Test pileup reweighter #
#############################

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
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from reweighting.pileupreweighter import PileupReweighter

# input arguments:
parser = argparse.ArgumentParser(description='Test pileup reweighter')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-s', '--sffile', required=True, type=os.path.abspath)
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

# make and evaluate a reweighter
reweighter = PileupReweighter(args.sffile, year)
weights = reweighter.weights(events)
weightsup = reweighter.weightsup(events)
weightsdown = reweighter.weightsdown(events)

# make a plot
fig, ax = plt.subplots()
xaxrange = (0.5, 2.)
nbins = 20
ax.hist(np.clip(weights, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.hist(np.clip(weightsup, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Up variation', color='cyan')
ax.hist(np.clip(weightsdown, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Down variation', color='fuchsia')
ax.set_xlabel('Per-event weight')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test.png')
