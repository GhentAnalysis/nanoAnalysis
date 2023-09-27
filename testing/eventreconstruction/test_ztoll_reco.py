#################################################
# Test Z boson reconstruction from lepton pairs #
#################################################

# imports
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(os.path.abspath('../../'))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from samples.sample import year_from_sample_name
from eventreconstruction.zreco import ZReco

# input arguments:
parser = argparse.ArgumentParser(description='Z boson reconstruction')
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
print('Performing lepton selection...')
muon_loose_mask = muonselection(events.Muon, selectionid='run2ul_loose')
electron_cleaning_mask = clean_electrons_from_muons(events.Electron, events.Muon[muon_loose_mask])
electron_loose_mask = ( 
  (electronselection(events.Electron, selectionid='run2ul_loose'))
  & electron_cleaning_mask )

# reconstruct Z bosons
print('Performing Z boson reconstruction...')
zreco = ZReco(events, halfwindow=10., samesign=False,
          electron_mask=electron_loose_mask, muon_mask=muon_loose_mask)
ee_candidates = zreco.get_ztoee_candidates()
mm_candidates = zreco.get_ztomm_candidates()
ll_candidates = zreco.get_ztoll_candidates()
best_mass = zreco.get_ztoll_best_mass()
has_candidate = zreco.has_ztoll_candidate()

# printouts for testing
doprint = False
if doprint:
    for i in range(100):
        print('---')
        print(ee_candidates[i])
        print(mm_candidates[i])
        print(ll_candidates[i])
        print(best_mass[i])
        print(has_candidate[i])

# make a plot
doplot = True
if doplot:
    mass = np.array(ak.fill_none(best_mass, -1))
    mass = mass[mass > 0.]
    fig,ax = plt.subplots()
    ax.hist(mass, range=(71,111), bins=30)
    fig.savefig('test.png')
