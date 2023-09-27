#############################################################
# Test b-tagging reweighter, specifically the normalization #
#############################################################

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
from reweighting.btagreweighter import BTagReweighter

# input arguments:
parser = argparse.ArgumentParser(description='Test b-tagging reweighter')
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

# calculate object masks
print('Performing jet selection...')
jet_mask = jetselection(events.Jet, selectionid='run2ul_default')
print('Number of events: {}'.format(len(events)))
print('Number of selected jets: {}'.format(ak.sum(jet_mask)))

# make a reweighter
reweighter = BTagReweighter(args.sffile)
systematics=['lf', 'hf']
print(reweighter)

# determine weights before normalization
weights = reweighter.weights(events, jet_mask=jet_mask)
print('Average weight before normalization: {}'.format(np.mean(weights)))
weightssys = {}
for sys in systematics:
    weightssys[sys] = {
      'up': reweighter.weightsup(events, jet_mask=jet_mask, unctype=sys),
      'down': reweighter.weightsdown(events, jet_mask=jet_mask, unctype=sys)
    }

# normalize reweighter
reweighter.set_normalization(events.Jet[jet_mask], unctypes=systematics)
print(reweighter)

# determine weights after normalization
reweighter.normalize = True
weights_norm = reweighter.weights(events, jet_mask=jet_mask)
print('Average weight after normalization: {}'.format(np.mean(weights_norm)))
weightssys_norm = {}
for sys in systematics:
    weightssys_norm[sys] = {
      'up': reweighter.weightsup(events, jet_mask=jet_mask, unctype=sys),
      'down': reweighter.weightsdown(events, jet_mask=jet_mask, unctype=sys)
    }

# make a plot
fig, axs = plt.subplots(ncols=2, figsize=(10,6))
xaxrange = (0., 3.)
nbins = 20
def plot(ax, weights, weightssys, text):
    ax.hist(np.clip(weights, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
    for sys in systematics:
        ax.hist(np.clip(weightssys[sys]['up'], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
          histtype='step', linewidth=2, color='r', alpha=0.2)
        ax.hist(np.clip(weightssys[sys]['down'], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
          histtype='step', linewidth=2, color='r', alpha=0.2)
    ax.set_xlabel('Per-event weight')
    ax.set_ylabel('Number of events')
    ax.grid()
    ax.legend()
    ax.text(0.05, 0.95, text, ha='left', va='top', transform=ax.transAxes)
plot(axs[0], weights, weightssys, 'Before normalization')
plot(axs[1], weights_norm, weightssys_norm, 'After normalization')
fig.savefig('test.png')
