##########################################
# Make cutflow using TTW event selection #
##########################################

# imports
import sys
import os
import argparse
from pathlib import Path
import ROOT
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
import tools.argparsetools as apt

# input arguments:
parser = argparse.ArgumentParser(description='Cutflow using TTW event selection')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('--skimmed', default=False, action='store_true')
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

# do event selection
print('Performing event selection...')
metfilter_mask = tst.pass_met_filters(events)
trigger_mask = tst.pass_any_lepton_trigger(events)
fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==2)
lowmass_mask = lst.pass_mll_lowmass_veto(events,
  electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
photon_mask = sst.pass_photon_overlap_removal(events, samplename=samplename)
tight_mask = lst.pass_tight_lepton_selection(events, 2, 'tight',
  electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
  electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.),
  electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
invmass = lst.get_lepton_invmass(events,
  electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
invmass_mask = (invmass > 30.)
ss_mask = (abs(ak.sum(ak.concatenate((events.Electron[electron_fo_mask].charge,
  events.Muon[muon_fo_mask].charge),axis=1),axis=1))==2)
zveto_mask = (
  (ak.sum(electron_fo_mask,axis=1)!=2)
  | (abs(invmass - 91.1876) > 10.) )
met_mask = ( events.MET.pt > 30. )
nbjet_mask = (ak.sum(bjet_loose_mask,axis=1) >= 2)
njet_mask = (ak.sum(jet_mask,axis=1) >= 3)

# aggregate masks
masks = {
  'MET filters': metfilter_mask,
  'trigger': trigger_mask,
  '2 FO leptons': fo_mask,
  'low mass veto': lowmass_mask,
  'photon overlap': photon_mask,
  '2 tight leptons': tight_mask,
  'pT thresholds': pt_mask,
  'invariant mass veto': invmass_mask,
  'same sign': ss_mask,
  'electron Z veto': zveto_mask,
  'MET': met_mask,
  'b-tagged jets': nbjet_mask,
  'jets': njet_mask
}

# printouts for testing
for i in range(100):
    print('---')
    print(events.Electron.pt[electron_fo_mask][i])
    print(events.Muon.pt[muon_fo_mask][i])
    print(pt_mask[i])

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
