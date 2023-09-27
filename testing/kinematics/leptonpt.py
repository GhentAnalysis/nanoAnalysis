##########################################
# Make cutflow using TTW event selection #
##########################################

# imports
import sys
import os
import argparse
import ROOT
import matplotlib.pyplot as plt
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.jetselection import jetselection
from objectselection.bjetselection import bjetselection
from objectselection.cleaning import clean_electrons_from_muons
from objectselection.cleaning import clean_jets_from_leptons
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
import eventselection.lepton_selection_tools as lst
import eventselection.sample_selection_tools as sst
import eventselection.trigger_selection_tools as tst
from plotting.singlehistplotter import plotsinglehistogram
import tools.argparsetools as apt

# input arguments:
parser = argparse.ArgumentParser(description='Make lepton pt distribution')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('--skimmed', default=False, action='store_true')
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make output directory
if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

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
muon_loose_mask = muonselection(events.Muon, selectionid='run2ul_loose')
muon_fo_mask = muonselection(events.Muon, selectionid='ttwloose_fo')
muon_tight_mask = muonselection(events.Muon, selectionid='ttwloose_tight')
electron_cleaning_mask = clean_electrons_from_muons(events.Electron, events.Muon[muon_loose_mask])
electron_loose_mask = ( 
  (electronselection(events.Electron, selectionid='run2ul_loose'))
  & electron_cleaning_mask )
electron_fo_mask = ( 
  (electronselection(events.Electron, selectionid='ttwloose_fo'))
  & electron_cleaning_mask )
electron_tight_mask = ( 
  (electronselection(events.Electron, selectionid='ttwloose_tight'))
  & electron_cleaning_mask )

# get lepton pt
leppts = {}
leppts['electron_all_pt'] = ak.flatten(events.Electron.pt)
leppts['electron_loose_pt'] = ak.flatten(events.Electron.pt[electron_loose_mask])
leppts['electron_fo_pt'] = ak.flatten(events.Electron.pt[electron_fo_mask])
leppts['electron_tight_pt'] = ak.flatten(events.Electron.pt[electron_tight_mask])
leppts['muon_all_pt'] = ak.flatten(events.Muon.pt)
leppts['muon_loose_pt'] = ak.flatten(events.Muon.pt[muon_loose_mask])
leppts['muon_fo_pt'] = ak.flatten(events.Muon.pt[muon_fo_mask])
leppts['muon_tight_pt'] = ak.flatten(events.Muon.pt[muon_tight_mask])

# make simple plots
for flavour in ['electron','muon']:
    for selectiontype in ['all','loose','fo','tight']:
        outputfile = os.path.join(args.outputdir,'{}_{}_pt.png'.format(flavour, selectiontype))
        hist = ROOT.TH1D('hist', 'hist', 20, 0., 150.)
        for el in leppts['{}_{}_pt'.format(flavour,selectiontype)]: hist.Fill(el)
        plotsinglehistogram(hist, outputfile,
                xaxtitle='Lepton p_{T} (GeV)', yaxtitle='Number of leptons',
                drawoptions='',
                extrainfos=[flavour, selectiontype])   
