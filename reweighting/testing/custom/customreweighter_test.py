##########################
# Test custom reweighter #
##########################

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
from reweighting.customreweighter import CustomReweighter


# define evaluator
def njets_evaluator(events, **kwargs):
    # initialize
    nominal = np.ones(len(events))
    up = np.ones(len(events))
    down = np.ones(len(events))
    # get jets
    if 'jet_mask' not in kwargs.keys():
        msg = 'ERROR: reweighter.weight was called without'
        msg += ' required argument "jet_mask".'
        raise Exception(msg)
    jets = events.Jet[kwargs['jet_mask']]
    njets = ak.num(jets) 
    # 20% uncertainty for 3 or more jets
    up = np.where(njets>=3, 1.2, up)
    down = np.where(njets>=3, 0.8, down)
    # 50% uncertainty for 6 or more jets
    up = np.where(njets>=6, 1.5, up)
    down = np.where(njets>=6, 0.5, down)
    return {'nominal': nominal, 'up': up, 'down': down}

    
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
print('Number of events: {}'.format(len(events)))
print('Number of selected jets: {}'.format(ak.sum(jet_mask)))

# make and evaluate a reweighter
reweighter = CustomReweighter(evaluator=njets_evaluator)
weights = reweighter.weights(events, jet_mask=jet_mask)
weightsup = reweighter.weightsup(events, jet_mask=jet_mask)
weightsdown = reweighter.weightsdown(events, jet_mask=jet_mask)

# calculate variable to plot
njets = ak.num(events.Jet[jet_mask])

# make a plot
fig, ax = plt.subplots()
xaxrange = (-0.5, 8.5)
nbins = int(xaxrange[1]-xaxrange[0])
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=weights,
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=weightsup,
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Up variation', color='cyan')
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=weightsdown,
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Down variation', color='fuchsia')
ax.set_xlabel('Number of jets')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test.png')
