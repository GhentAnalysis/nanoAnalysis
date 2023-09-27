###############################
# Test lepton cone correction #
###############################

# imports
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.conecorrectionfactors import electron_cone_correction_factor
from objectselection.conecorrectionfactors import muon_cone_correction_factor
from preprocessing.preprocessor import PreProcessor
from preprocessing.conecorrection import lepton_cone_correction
from samples.sample import year_from_sample_name

# input arguments:
parser = argparse.ArgumentParser(description='Test lepton cone correction')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
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
electron_loose_mask = electronselection(events.Electron, selectionid='run2ul_loose')
electron_fo_mask = electronselection(events.Electron, selectionid='ttwloose_fo')
electron_tight_mask = electronselection(events.Electron, selectionid='ttwloose_tight')

# calculate additional object masks
muon_fonottight_mask = ((muon_fo_mask) & (~muon_tight_mask))
electron_fonottight_mask = ((electron_fo_mask) & (~electron_tight_mask))

# calculate event masks
event_fomuon_mask = (ak.sum(muon_fonottight_mask, axis=1)>=1)
event_tightmuon_mask = (ak.sum(muon_tight_mask, axis=1)>=1)
event_foelectron_mask = (ak.sum(electron_fonottight_mask, axis=1)>=1)
event_tightelectron_mask = (ak.sum(electron_tight_mask, axis=1)>=1)

# pt before correction
print('=== Before correction ===')
print('FO not tight muons:')
print(events.Muon[muon_fonottight_mask][event_fomuon_mask].pt)
print('tight muons:')
print(events.Muon[muon_tight_mask][event_tightmuon_mask].pt)
print('FO not tight electrons:')
print(events.Electron[electron_fonottight_mask][event_foelectron_mask].pt)
print('tight electrons:')
print(events.Electron[electron_tight_mask][event_tightelectron_mask].pt)

# do correction
lepton_cone_correction(events,
  electron_fo_mask=electron_fo_mask, electron_tight_mask=electron_tight_mask,
  muon_fo_mask=muon_fo_mask, muon_tight_mask=muon_tight_mask,
  electron_correctionfactor=electron_cone_correction_factor('ttwloose'),
  muon_correctionfactor=muon_cone_correction_factor('ttwloose')
)

# pt after correction
print('=== After correction ===')
print('FO not tight muons:')
print(events.Muon[muon_fonottight_mask][event_fomuon_mask].pt)
print(events.Muon[muon_fonottight_mask][event_fomuon_mask].ptNoCone)
print('tight muons:')
print(events.Muon[muon_tight_mask][event_tightmuon_mask].pt)
print(events.Muon[muon_tight_mask][event_tightmuon_mask].ptNoCone)
print('FO not tight electrons:')
print(events.Electron[electron_fonottight_mask][event_foelectron_mask].pt)
print(events.Electron[electron_fonottight_mask][event_foelectron_mask].ptNoCone)
print('tight electrons:')
print(events.Electron[electron_tight_mask][event_tightelectron_mask].pt)
print(events.Electron[electron_tight_mask][event_tightelectron_mask].ptNoCone)

# try another cone correction (should raise exception)
lepton_cone_correction(events,
    electron_fo_mask=electron_fo_mask, electron_tight_mask=electron_tight_mask,
    muon_fo_mask=muon_fo_mask, muon_tight_mask=muon_tight_mask)
