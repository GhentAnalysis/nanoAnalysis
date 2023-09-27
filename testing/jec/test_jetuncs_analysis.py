#####################################################
# Test JEC uncertainties and related jet operations #
#####################################################

# imports
import sys
import os
import argparse
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.jetselection import jetselection
from objectselection.jetuncs import get_varied_jets, get_varied_met
from objectselection.jetuncs import get_available_jec_variations
from objectselection.jetuncs import get_available_jer_variations
from samples.sample import year_from_sample_name


def pass_eventselection(events, jet_mask=None):
    ### dummy event selection
    njet_mask = (ak.sum(jet_mask,axis=1) >= 5)
    # aggregate masks
    masks = {
      'Jets': njet_mask
    }
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def calculate_event_variables(events, jet_mask=None):
    ### calculate dummy event variables
    res = {}
    # initializations
    nevents = ak.count(events.event)
    # total yield (fixed arbitrary value)
    res['yield'] = np.ones(nevents)*0.5
    # get object collections
    jets = events.Jet[jet_mask]
    # number of jets
    njets = ak.sum(jet_mask, axis=1)
    res['nJets'] = njets
    # jet pt and eta
    jet_pt = np.array(ak.fill_none(ak.pad_none(jets.pt, 2, axis=1, clip=True), 0.))
    res['jetPtLeading'] = jet_pt[:,0]
    return res


# input arguments:
parser = argparse.ArgumentParser(description='Test JEC uncertainties')
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

# print available variations
all_jecs = get_available_jec_variations(events)
print(all_jecs)
all_jers = get_available_jer_variations(events)
print(all_jers)

# initializations
sources = all_jecs + all_jers + ['unclustEn']
variations = ['nominal']
for source in sources:
    variations.append(source+'Up')
    variations.append(source+'Down')
res = {}

# store original jets
jets_nominal = events.Jet

# test approach 1: replace events.Jet collection
# advantage: easy syntax, one-time calculation, no need for additional arguments
# disadvantage: messing with events object might introduce nasty bugs

# loop over variations
for variation in variations:
    print('Now running on variation {}'.format(variation))

    # get varied jet collection
    jets = get_varied_jets(events, variation)
    events.Jet = jets
    events['Jet'] = jets

    # do object and event selection
    jet_mask = jetselection(jets, selectionid='run2ul_default')
    events_mask = pass_eventselection(events, jet_mask=jet_mask)

    # calculate event variables
    variables = calculate_event_variables(events, jet_mask=jet_mask)
    variables = {key: val[events_mask] for key,val in variables.items()}
    res[variation] = variables

    # printouts for comparison
    print('Selected {} events'.format(ak.sum(events_mask)))
    print('Average number of jets: {}'.format(np.mean(variables['nJets'])))
    print('Average leading jet pt: {}'.format(np.mean(variables['jetPtLeading'])))

    # re-set original jets
    events.Jet = jets_nominal
    events['Jet'] = jets_nominal

# make plots
for source in sources:
    fig, ax = plt.subplots()
    xaxrange = (3.5, 9.5)
    plotvar = 'nJets'
    nbins = int(xaxrange[1]-xaxrange[0])
    ax.hist(np.clip(res['nominal'][plotvar], xaxrange[0], xaxrange[1]),
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Nominal', color='dodgerblue')
    ax.hist(np.clip(res[source+'Up'][plotvar], xaxrange[0], xaxrange[1]),
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Up variation', color='cyan')
    ax.hist(np.clip(res[source+'Down'][plotvar], xaxrange[0], xaxrange[1]),
      bins=nbins, range=xaxrange,
      histtype='step', linewidth=2, label='Down variation', color='fuchsia')
    ax.set_xlabel('Number of jets')
    ax.set_ylabel('Number of events')
    ax.grid()
    ax.legend()
    fig.savefig('test_{}.png'.format(source))
