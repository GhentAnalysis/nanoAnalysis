####################################
# Test script for object selection #
####################################
# in particular, for synchronization of lepton ID with miniAOD based framework


# import python modules
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
# import framework modules
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from samples.sample import dtype_from_sample_name
import eventselection.lepton_selection_tools as lst
import tools.argparsetools as apt


if __name__=='__main__':


  # input arguments:
  parser = argparse.ArgumentParser(description='Synchronization of lepton ID')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-n', '--nentries', type=int, default=-1)
  parser.add_argument('--skimmed', default=False, action='store_true')
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))
  sys.stdout.flush()

  # make NanoEvents array
  print('Loading events from input file...')
  year = year_from_sample_name(args.inputfile)
  dtype = dtype_from_sample_name(args.inputfile)
  samplename = os.path.basename(args.inputfile)
  events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year, 'samplename': samplename, 'dtype': dtype}
  ).events()
  nevents = ak.count(events.event)
  print('Number of events in input file: {}'.format(nevents))
  sys.stdout.flush()

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

  # print lepton variables
  for i in range(nevents):
    print('-----')
    print('Event: {}'.format(events.event[i]))
    print('Electrons:')
    print(events.Electron.pt[i])
    print(events.Electron.eta[i])
    print(events.Electron.jetBTagDeepFlavor[i])
    print(events.Electron.jetPtRatio[i])
    print(events.Electron.mvaFall17V2noIso_WPL[i])
    print('Muons:')
    print(events.Muon.pt[i])
    print(events.Muon.eta[i])
    print(events.Muon.jetBTagDeepFlavor[i])
    print(events.Muon.jetPtRatio[i])

  # calculate object masks
  print('Performing lepton selection...')
  sys.stdout.flush()
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
