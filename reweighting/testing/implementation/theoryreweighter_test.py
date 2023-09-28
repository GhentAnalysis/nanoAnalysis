##########################
# Test theory reweighter #
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
from samples.sampleweights import SampleWeights
from samples.sample import year_from_sample_name
from reweighting.implementation.run2ulreweighter import get_run2ul_reweighter
from objectselection.jetselection import jetselection


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
jet_mask = jetselection(events.Jet, selectionid='run2ul_default')

# make theory reweighter
sampleweights = SampleWeights(args.inputfile)
reweighter = get_run2ul_reweighter(year, sampleweights)

# disable some uncertainties
keep_reweighters = ['pdfacceptance']
#keep_reweighters = ['pdfnorm']
#keep_reweighters = ['scaleacceptance']
#keep_reweighters = ['scalenorm']
#keep_reweighters = ['ps']
all_names = list(reweighter.reweighters.keys())
for name in all_names:
    if name not in keep_reweighters: reweighter.remove_reweighter(name)

# determine args and uncertainties
args = reweighter.get_args()
uncs = reweighter.get_unctypes()
variations = reweighter.get_variations()

# printouts for testing
print('All arguments:')
print(args)
print('All uncertainties:')
print(uncs)
print('All variations')
print(variations)

# get all weights
print('Calculating weights')
res = {}
res['nominal'] = reweighter.weights(events)
for name in reweighter.reweighters.keys():
    for variation in variations[name]:
        key = '_'.join([name, str(variation)])
        res[key] = reweighter.weightsvar(events, name, variation)

# get variable to plot
njets = ak.num(events.Jet[jet_mask])

# make a plot of the weights
print('Making weight plot')
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
print('Making variable plot')
fig, ax = plt.subplots()
xaxrange = (-0.5, 8.5)
nbins = int(xaxrange[1]-xaxrange[0])
for sys in res.keys():
    if sys=='nominal': continue
    ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=res[sys],
          bins=nbins, range=xaxrange,
          histtype='step', linewidth=2, color='r', alpha=0.2)
ax.hist(np.clip(njets, xaxrange[0], xaxrange[1]), weights=res['nominal'],
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
ax.set_xlabel('Number of jets')
ax.set_ylabel('Number of events')
ax.grid()
ax.legend()
fig.savefig('test2.png')
