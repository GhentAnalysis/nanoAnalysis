################################################
# Script for making ROC curves for lepton MVAs #
################################################
# part 1: filling histograms and writing to file

# imports
import sys
import os
from pathlib import Path
import argparse
import numpy as np
import uproot
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# local imports
sys.path.append(str(Path(__file__).parents[2]))
from preprocessing.preprocessor import PreProcessor
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from samples.sample import year_from_sample_name


# input arguments:
parser = argparse.ArgumentParser(description='Make ROC curves for lepton MVAs')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# define settings for which to make ROC curves
#leptonid_pre = 'topmvav2_loose'
leptonid_pre = 'run2ul_loose'
#leptonid_post = 'run2ul_loose'
leptonid_post = None
nbins = 5000
binrange = (-1.,1.)
flavours = ['electron', 'muon']
mvavars = ['mvaTTH', 'mvaTOP']
ptcuts = [(10.,0), (25.,0), (10.,25.)]

# make NanoEvents array
year = year_from_sample_name(args.inputfile)
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year}
).events()

# calculate additional needed variables
if( not hasattr(events.Electron, 'mvaTOP')
    or not hasattr(events.Muon, 'mvaTOP') ):
    preprocessor = PreProcessor()
    preprocessor.process(events,
      leptongenvariables=['isPrompt'],
      leptonvariables=['jetPtRatio','jetBTagDeepFlavor'],
      topmvavariable='mvaTOP', topmvaversion='ULv1')

# select leptons with appropriate pre-selection
# (leptons not selected as this stage are not taken into account)
selected_muons = events.Muon[muonselection(events.Muon, selectionid=leptonid_pre)]
selected_electrons = events.Electron[ 
  ( (electronselection(events.Electron, selectionid=leptonid_pre))
    & (clean_electrons_from_muons(events.Electron, selected_muons)) )
]

# put mva score of leptons not passing additional selections to minimum value
if leptonid_post is not None:
    electron_post_mask = electronselection(selected_electrons, selectionid=leptonid_post)
    muon_post_mask = muonselection(selected_muons, selectionid=leptonid_post)
    selected_electrons.mvaTTH = ak.where(electron_post_mask, selected_electrons.mvaTTH, -1)
    if 'mvaTOP' in mvavars:
        selected_electrons.mvaTOP = ak.where(electron_post_mask, selected_electrons.mvaTOP, -1)
    selected_muons.mvaTTH = ak.where(muon_post_mask, selected_muons.mvaTTH, -1)
    if 'mvaTOP' in mvavars:
        selected_muons.mvaTOP = ak.where(muon_post_mask, selected_muons.mvaTOP, -1)
    # implementation note: how to generalize the above to an arbitrary list of mvavars?
    # cannot use getattr, since it throws the error 'cannot assign a function call'

# define info dict
mvascores = {}

# loop over flavours
for flavour in flavours:
    # do flavour selection
    if flavour=='electron': selected_leptons = selected_electrons
    elif flavour=='muon': selected_leptons = selected_muons
    else: raise Exception('ERROR: flavour {} not recognized.'.format(flavour))
    # loop over pt cuts
    for ptcut in ptcuts:
        # define pt mask
        ptmask = (selected_leptons.pt>0)
        if ptcut[0]>0: ptmask = ptmask & (selected_leptons.pt>ptcut[0])
        if ptcut[1]>0: ptmask = ptmask & (selected_leptons.pt<ptcut[1])
        # define pt name
        ptname = 'pt{}to{}'.format(int(ptcut[0]), int(ptcut[1]) if ptcut[1]>0 else 'Inf')
        # loop over mvas
        for mvavar in mvavars:
            # loop over prompt/nonprompt
            for p in ['prompt','nonprompt']:
                # define prompt/nonprompt mask
                pmask = (selected_leptons.isPrompt)
                if p=='nonprompt': pmask = ~pmask
                # define a name for this instance
                name = '_'.join([flavour,ptname,mvavar,p])
                # get and bin the values
                values = getattr(selected_leptons[ptmask & pmask], mvavar)
                values = np.array(ak.flatten(values))
                hist = np.histogram(values, bins=nbins, range=binrange)
                # add to info dict
                mvascores[name] = hist

# write to output file
with uproot.recreate(args.outputfile) as f:
    for key,val in mvascores.items():
        f[key] = val
