###########################
# Test prefire reweighter #
###########################

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
from reweighting.prefirereweighter import PrefireReweighter

# input arguments:
parser = argparse.ArgumentParser(description='Test prefire reweighter')
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

# make and evaluate a reweighter
ptypes = ['all', 'ecal', 'muon']
data = {}
for ptype in ptypes:
    unctype = None
    if ptype=='muon': unctype='Syst'
    reweighter = PrefireReweighter(prefiretype=ptype)
    weights = reweighter.weights(events)
    weightsup = reweighter.weightsup(events, unctype=unctype)
    weightsdown = reweighter.weightsdown(events, unctype=unctype)
    data[ptype] = {'nominal': weights, 'up': weightsup, 'down': weightsdown}

# make a plot
fig, axs = plt.subplots(ncols=3, figsize=(15,6))
xaxrange = (0.97, 1.01)
nbins = 20
def plot(ax, weights, weightsup, weightsdown, text):
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
    ax.text(0.05, 0.95, text, ha='left', va='top', transform=ax.transAxes)
for i, ptype in enumerate(ptypes):
    plot(axs[i], data[ptype]['nominal'], data[ptype]['up'], data[ptype]['down'],
         'Prefire type: {}'.format(ptype))
fig.savefig('test.png')
