###############################
# Test electron ID reweighter #
###############################

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
from objectselection.electronselection import electronselection
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from reweighting.electronidreweighter import ElectronIDReweighter

# input arguments:
parser = argparse.ArgumentParser(description='Test electron ID reweighter')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-s', '--sffile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('--skimmed', default=False, action='store_true')
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

# calculate additional variables
if args.skimmed:
    preprocessor = PreProcessor()
    preprocessor.process(events)
else:
    preprocessor = PreProcessor()
    preprocessor.process(events,
        leptongenvariables=[
          'isPrompt',
        ],
        leptonvariables=[
          'jetPtRatio',
          'jetBTagDeepFlavor'
        ],
        topmvavariable='mvaTOP', topmvaversion='ULv1',
        dotriggers=True
    )

# calculate object masks
print('Performing lepton selection...')
electron_mask = electronselection(events.Electron, selectionid='ttwloose_tight')
print('Number of events: {}'.format(len(events)))
print('Number of selected leptons: {}'.format(ak.sum(electron_mask)))

# make and evaluate a reweighter
reweighter = ElectronIDReweighter(args.sffile)
weights = reweighter.weights(events, electron_mask=electron_mask)
weightsup = reweighter.weightsup(events, electron_mask=electron_mask)
weightsdown = reweighter.weightsdown(events, electron_mask=electron_mask)

# printouts for testing
print(weights)

# make a plot
fig, ax = plt.subplots()
xaxrange = (0.8, 1.2)
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
