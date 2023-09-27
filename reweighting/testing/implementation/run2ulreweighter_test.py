############################
# Test combined reweighter #
############################

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
from samples.sample import year_from_sample_name
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
reweighter = get_run2ul_reweighter(year)
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
uncs = reweighter.get_uncertainties()

# printouts for testing
print('All arguments:')
print(args)
print('All uncertainties:')
print(uncs)

# get all weights
res = {}
res['nominal'] = reweighter.weights(events, **kwargs)
for name in reweighter.reweighters.keys():
    unctypes = uncs[name]
    if unctypes is None:
        up = reweighter.weightsup(events, name, **kwargs)
        down = reweighter.weightsdown(events, name, **kwargs)
        res[name] = (up, down)
    else:
      for unctype in unctypes:
        key = '_'.join([name, unctype])
        up = reweighter.weightsup(events, name, unctype=unctype, **kwargs)
        down = reweighter.weightsdown(events, name, unctype=unctype, **kwargs)
        res[key] = (up, down)

# get variable to plot
njets = ak.num(events.Jet[jet_mask])

# make a plot of the weights
fig, ax = plt.subplots()
xaxrange = (0., 3.)
nbins = 20
for sys in res.keys():
    if sys=='nominal': continue
    for i in [0,1]:
        ax.hist(np.clip(res[sys][i], xaxrange[0], xaxrange[1]), bins=nbins, range=xaxrange,
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
    for i in [0,1]:
        ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=res[sys][i],
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
