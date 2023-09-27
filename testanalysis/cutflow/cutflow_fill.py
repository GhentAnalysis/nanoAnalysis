###################################################################
# Perform TTW event selection and calculate interesting variables #
###################################################################
# input: a single nanoAOD file.
# output: a file containing a ROOT tree with event variables


# import python modules
import sys
import os
import argparse
from pathlib import Path
import ROOT
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
# import framework modules
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.jetselection import jetselection
from objectselection.bjetselection import bjetselection
from objectselection.cleaning import clean_electrons_from_muons
from objectselection.cleaning import clean_jets_from_leptons
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from samples.sample import dtype_from_sample_name
from samples.sampleweights import SampleWeights
import eventselection.lepton_selection_tools as lst
import eventselection.sample_selection_tools as sst
import eventselection.trigger_selection_tools as tst
import tools.argparsetools as apt
from tools.readfakeratetools import readfrmapfromfile
from tools.readchargefliptools import readcfmapfromfile
# import local modules
sys.path.append(os.path.abspath('../eventselections'))
from eventselections import pass_event_selection


if __name__=='__main__':

  # input arguments:
  parser = argparse.ArgumentParser(description='Perform cutflow analysis')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
  parser.add_argument('-n', '--nentries', type=int, default=-1)
  parser.add_argument('-s', '--eventselection', required=True)
  parser.add_argument('-t', '--selectiontype', default='tight')
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

  # calculate object masks
  print('Performing object selection...')
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
  leptonsforcleaningjets = ak.with_name(ak.concatenate(
    (events.Electron[electron_fo_mask],events.Muon[muon_fo_mask]), axis=1),
    'PtEtaPhiMCandidate')
  jet_cleaning_mask = clean_jets_from_leptons(events.Jet, leptonsforcleaningjets)
  jet_mask = (
    jetselection(events.Jet, selectionid='run2ul_default')
    & jet_cleaning_mask )
  bjet_loose_mask = (
    jet_mask 
    & bjetselection(events.Jet, year=year, algo='deepflavor', level='loose') )

  # loop over event selections and selection types
  print('Performing event selection for {} / {}...'.format(args.eventselection, args.selectiontype))
  sys.stdout.flush()

  # do event selection
  masks = pass_event_selection(events,
            args.eventselection,
            selectiontype=args.selectiontype,
            electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
            electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask,
            jet_mask=jet_mask, bjet_mask=bjet_loose_mask,
            cutflow=True)

  # printouts for testing
  '''concatfo = ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1)
  sumfo = ak.sum(concatfo,axis=1)
  maskfo = (sumfo==2)
  concatt = ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1)
  sumt = ak.sum(concatt,axis=1)
  maskt = (sumt<2)
  for i in range(20):
    print('-----')
    print(events.event[i])
    #print(concatfo[i])
    #print(sumfo[i])
    #print(maskfo[i])
    print(concatt[i])
    print(sumt[i])
    print(maskt[i])'''

  # more printouts for testing
  '''maskkeys = list(masks.keys())
  totalmask = (events.event > 0)
  #with open('temp_test_Nomask.txt', 'w') as f:
  #  for i in range(len(events)):
  #    f.write('{}\n'.format(events.event[i]))
  for i,key in enumerate(maskkeys):
    totalmask = totalmask & masks[key]
    if i<5: continue
    tempevents = events[totalmask]
    with open('temp_test_{}.txt'.format(key.replace(' ','')), 'w') as f:
      for i in range(len(tempevents)):
        f.write('{} {} {}\n'.format(tempevents.run[i], tempevents.luminosityBlock[i], tempevents.event[i]))'''

  # do printouts
  nmasks = len(masks)
  totalmask = (events.event >= 0)
  totalremaining = ak.sum(totalmask)
  failatstep = []
  for key, val in masks.items():
    totalmask = totalmask & val
    totalremaining_new = ak.sum(totalmask)
    failatstep.append(totalremaining-totalremaining_new)
    totalremaining = totalremaining_new
    print('{}: {} events remaining'.format(key, totalremaining))

  # make a histogram for plotting
  # implementation note: the choice of histogram is simply to be
  # fully consistent with earlier nanoAOD tests with ewkino framework,
  # for easier synchronization checks (and to use the same plotting function).
  cutflowhist = ROOT.TH1D('cutflowhist', 'cutflowhist;cut;events', len(masks)+1, -0.5, len(masks)+0.5)
  for i, (key,val) in enumerate(masks.items()):
    cutflowhist.GetXaxis().SetBinLabel(i+1, 'Fail '+key)
    cutflowhist.SetBinContent(i+1, failatstep[i])
  cutflowhist.GetXaxis().SetBinLabel(len(masks)+1, 'Pass')
  cutflowhist.SetBinContent(len(masks)+1, totalremaining)
  f = ROOT.TFile.Open(args.outputfile, 'recreate')
  cutflowhist.Write()
  f.Close()
