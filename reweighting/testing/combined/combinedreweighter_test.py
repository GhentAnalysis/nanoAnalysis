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
from objectselection.jetselection import jetselection
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from reweighting.muonidreweighter import MuonIDReweighter
from reweighting.electronidreweighter import ElectronIDReweighter
from reweighting.pileupreweighter import PileupReweighter
from reweighting.prefirereweighter import PrefireReweighter
from reweighting.combinedreweighter import CombinedReweighter


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

# make combined reweighter
reweighter = CombinedReweighter()
wfile = os.path.join('../../data/pileup', 'puWeights_{}.json'.format(year))
reweighter.add_reweighter('pileup', PileupReweighter(wfile, year))
reweighter.add_reweighter('prefire', PrefireReweighter(prefiretype='all'))

# printouts for testing
print(reweighter.get_unctypes())
print(reweighter.get_variations())

# calculate total weights
res1 = reweighter.allweights(events, wtype='total', verbose=True)
print(res1.keys())
print(res1['pileup_up'])

# calculate individual weights
res2 = reweighter.allweights(events, wtype='individual', verbose=True)
print(res2.keys())
print(res2['pileup_up'])
tot = np.multiply(np.divide(res2['nominal'],res2['pileup']),res2['pileup_up'])
print(tot)
