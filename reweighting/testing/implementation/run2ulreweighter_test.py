############################
# Test combined reweighter #
############################

# imports
import sys
import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[3]))
from samples.sample import year_from_sample_name
from samples.sampleweights import SampleWeights
from reweighting.implementation.run2ulreweighter import get_run2ul_reweighter
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.jetselection import jetselection
from objectselection.bjetselection import bjetselection


# input arguments:
parser = argparse.ArgumentParser(description='Test combined reweighter')
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
print('Performing object selection...')
muon_mask = muonselection(events.Muon, selectionid='run2ul_loose')
electron_mask = electronselection(events.Electron, selectionid='run2ul_loose')
jet_mask = jetselection(events.Jet, selectionid='run2ul_default')
bjet_mask = (jet_mask
    & bjetselection(events.Jet, year=year, algo='deepflavor', level='loose') )

# make combined reweighter
sampleweights = SampleWeights(args.inputfile)
reweighter = get_run2ul_reweighter(year, sampleweights)
kwargs = ({
  'electron_mask': electron_mask,
  'muon_mask': muon_mask,
  'jet_mask': jet_mask,
  'bjet_mask': bjet_mask
})

# remove some reweighters that are not relevant
reweighter.remove_reweighter('njets')
reweighter.remove_reweighter('nbjets')

# normalize b-tag reweighter
btagreweighter = reweighter.reweighters['btagging']
btagreweighter.set_normalization(events.Jet[jet_mask], unctypes='all')

# determine args and uncertainties
args = reweighter.get_args()
uncs = reweighter.get_unctypes()
variations = reweighter.get_variations()

# printouts for testing
print('All arguments:')
print(args)
print('All uncertainties:')
print(uncs)
print('All variations:')
print(variations)

# get all weights
res = {}
res['nominal'] = reweighter.weights(events, **kwargs)
start_time = time.time()
approach = 3
for name in reweighter.reweighters.keys():
    print('now calculating {} variations...'.format(name))
    sys.stdout.flush()
    if approach==1:
        # approach 1: use weightsvar
        for variation in variations[name]:
            key = '_'.join([name, str(variation)])
            res[key] = reweighter.weightsvar(events, name, variation, **kwargs)
    if approach==2:
        # approach 2: use singleweightsvar
        # (much faster!)
        weights_single_nominal = reweighter.singleweights(events, name, **kwargs)
        zero_inds = np.nonzero(weights_single_nominal==0)
        weights_single_nominal[zero_inds] = 1.
        weights_withoutsingle = np.divide(res['nominal'], weights_single_nominal)
        for variation in variations[name]:
            key = '_'.join([name, str(variation)])
            weights_single_variation = reweighter.singleweightsvar(events, name, variation, **kwargs)
            res[key] = np.multiply(weights_withoutsingle, weights_single_variation)
if approach==3:
    # approach 3: same as above but wrapped in a convenience function in CombinedReweighter
    res = reweighter.allweights(events, wtype='total', verbose=True, **kwargs)
finish_time = time.time()
print('Calculation of systematics took {:.2f} seconds'.format(finish_time - start_time))

# get variable to plot
njets = ak.num(events.Jet[jet_mask])

# make a plot of the weights
fig, ax = plt.subplots()
xaxrange = (0., 3.)
nbins = 20
for sys in res.keys():
    if sys=='nominal': continue
    ax.hist(np.clip(res[sys], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
          histtype='step', linewidth=2, color='r', alpha=0.2)
ax.hist(np.clip(res['nominal'], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
  histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.set_xlabel('Per-event weight')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test1.png')

# make a plot of a variable
fig, ax = plt.subplots()
xaxrange = (-0.5, 8.5)
nbins = int(xaxrange[1]-xaxrange[0])
for sys in res.keys():
    if sys=='nominal': continue
    ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=res[sys],
          bins=nbins, range=xaxrange,
          histtype='step', linewidth=2, color='r', alpha=0.2)
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=np.ones(len(njets)),
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Unweighted', color='b')
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=res['nominal'],
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.set_xlabel('Number of jets')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test2.png')
