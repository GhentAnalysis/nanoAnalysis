##################################################################
# Same as btagreweighter_test.py but using a combined reweighter #
##################################################################
# The goal is to check that the output is exactly the same in both cases.

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
from reweighting.combinedreweighter import CombinedReweighter
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

# make and evaluate a reweighter
reweighter = CombinedReweighter()
reweighter.add_reweighter('btagging', BTagReweighter(args.sffile))
weights = reweighter.weights(events, jet_mask=jet_mask)
weightsup = reweighter.weightsup(events, 'btagging', jet_mask=jet_mask)
weightsdown = reweighter.weightsdown(events, 'btagging', jet_mask=jet_mask)
weightssplit = {}
systematics = reweighter.get_reweighter('btagging').get_unctypes()
for sys in systematics:
    weightssplit[sys] = {
      'up': reweighter.weightsup(events, 'btagging', jet_mask=jet_mask, unctype=sys),
      'down': reweighter.weightsdown(events, 'btagging', jet_mask=jet_mask, unctype=sys)
    }

# make a plot
fig, ax = plt.subplots()
xaxrange = (0., 3.)
nbins = 20
ax.hist(np.clip(weights, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.hist(np.clip(weightsup, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Up variation', color='cyan')
ax.hist(np.clip(weightsdown, xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Down variation', color='fuchsia')
for sys in systematics:
    ax.hist(np.clip(weightssplit[sys]['up'], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, color='r', alpha=0.2)
    ax.hist(np.clip(weightssplit[sys]['down'], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, color='r', alpha=0.2)
ax.set_xlabel('Per-event weight')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test2.png')
