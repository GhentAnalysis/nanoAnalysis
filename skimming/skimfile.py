##############################
# Skim a single nanoAOD file #
##############################
# The input is supposed to be a nanoAOD ROOT file.
# The output is a file of approximately the same format,
# with the following modifications (all optional):
# - addition of branches (e.g. TOP lepton MVA score for electrons and muons)
# - removal of unused branches (e.g. FatJet)
# - event selection

# imports
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[1]))
from skimming.nanoeventswriter import NanoEventsWriter, AuxiliaryTreeWriter
from skimming.skimselection import skimselection
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from samples.sample import dtype_from_sample_name
import tools.argparsetools as apt

# print starting tag
sys.stderr.write('###starting###\n')

# input arguments:
parser = argparse.ArgumentParser(description='Skim a single nanoAOD file')
parser.add_argument('-i', '--inputfile', required=True, type=apt.path_or_das)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('-l', '--leptonselection', default=None)
parser.add_argument('-s', '--skimselection', default=None)
parser.add_argument('-d', '--dropbranches', default=None)
parser.add_argument('--selectfirst', default=False, action='store_true',
                    help='Do selection first, before calculating additional variables;'
                        +' this might speed up the skimming,'
                        +' but will not work if additional variables are needed for selection.'
                        +' WARNING: do not use yet, gives errors on writing...')
parser.add_argument('--twosteps', default=False, action='store_true',
                    help='Do selection and additional variables as two separate steps,'
                        +' with an intermediate temporary file.'
                        +' The effect is similar as --selectfirst, but with intermediate writing.')
parser.add_argument('--compressionlevel', default=9, type=int,
                    help='Strength of compression (0 = no compression, 9 = maximal compression).')
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
print('Loading events from input file...')
year = year_from_sample_name(args.inputfile)
dtype = dtype_from_sample_name(args.inputfile)
print('Sample is found to be {} {}'.format(year,dtype))
uproot.open.defaults["xrootd_handler"] = uproot.MultithreadedXRootDSource
uproot.open.defaults["timeout"] = 360
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year, 'dtype': dtype}
).events()
print('Number of events in input file: {}'.format(ak.count(events.event)))

def preprocess(nanoevents):
    # calculate additional variables
    preprocessor = PreProcessor()
    leptongenvariables = [
        'isPrompt',
        'matchPdgId',
        'isChargeFlip'
    ]
    if nanoevents.metadata['dtype']=='data':
        leptongenvariables = None
    preprocessor.process(nanoevents,
      leptongenvariables=leptongenvariables,
      leptonvariables=[
        'jetPtRatio',
        'jetBTagDeepFlavor'
      ],
      topmvavariable='mvaTOP', topmvaversion='ULv1',
      dotriggers=True
    )

if not (args.selectfirst or args.twosteps):
    # calculate additonal variables
    print('Calculating additional variables on top of nanoAOD...')
    preprocess(events)

# calculate object masks
print('Performing lepton selection...')
muon_mask = muonselection(events.Muon, selectionid=args.leptonselection)
electron_mask = ( (electronselection(events.Electron, selectionid=args.leptonselection))
                  & (clean_electrons_from_muons(events.Electron, events.Muon[muon_mask])) )

# do event selection
print('Performing event selection...')
selected_events = events[
  skimselection(events, selectionid=args.skimselection, 
                muon_mask=muon_mask, electron_mask=electron_mask)
]
nselected_events = ak.count(selected_events.event)
print('Number of events after skim selection: {}'.format(nselected_events))

# throw error if number of selected events is zero
# (an empty tree apparently gives errors in writing by uproot;
#  to solve later but assume for now this will not happen significantly)
if nselected_events==0:
    msg = 'ERROR: number of selected events is zero, cannot write tree.'
    raise Exception(msg)

if args.selectfirst:
    # calculate additional variables
    print('Calculating additional variables on top of nanoAOD...')
    preprocess(selected_events)

if args.twosteps:
    # write events to temporary file
    # note: set compression to 0 for speed, since the file is temporary anyway.
    print('Writing events to temporary file...')
    tempfile = args.outputfile.replace('.root','_temp.root')
    writer = NanoEventsWriter()
    writer.write( selected_events, tempfile,
      compressionlevel=0, drop=args.dropbranches )
    # read events from temporary file
    print('Loading events from temporary file...')
    selected_events = NanoEventsFactory.from_root(
        tempfile,
        entry_stop=args.nentries if args.nentries>=0 else None,
        schemaclass=NanoAODSchema,
        metadata={'year': year, 'dtype': dtype}
    ).events()
    print('Number of events in input file: {}'.format(ak.count(selected_events.event)))
    # calculate additional variables
    print('Calculating additional variables on top of nanoAOD...')
    preprocess(selected_events)

# write to a file
print('Writing events to output file...')
writer = NanoEventsWriter()
writer.write( selected_events, args.outputfile,
  compressionlevel=args.compressionlevel, drop=args.dropbranches )

# copy auxiliary trees to new file
print('Writing auxiliary trees to output file...')
auxwriter = AuxiliaryTreeWriter()
auxwriter.write(args.inputfile, args.outputfile)

if args.twosteps:
    # delete temporary file
    print('Deleting temporary file...')
    os.system('rm {}'.format(tempfile))

# print done tag
sys.stderr.write('###done###\n')
