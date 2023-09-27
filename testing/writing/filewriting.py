######################################################################
# Basic test script for writing a NanoEvents array to a NanoAOD file #
######################################################################

# imports
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
sys.path.append(str(Path(__file__).parents[2]))
from skimming.nanoeventswriter import NanoEventsWriter
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.cleaning import clean_electrons_from_muons
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name

# input arguments:
parser = argparse.ArgumentParser(description='Test file writing and re-reading')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('-m', '--mode', choices=['write','read'], default='write')
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make NanoEvents array
year = None
if args.mode=='write': year = year_from_sample_name(args.inputfile)
events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year}
).events()
print('Basic properties of "events" array:')
print(events)
print(type(events))
print('Number of events: {}'.format(ak.count(events.event)))

if args.mode=='read':
    # print out all the fields
    print('Contents of "events" array:')
    for field in sorted(events.fields):
        print('  {}'.format(field))
        for subfield in sorted(events[field].fields):
            print('    {}'.format(subfield))
    sys.exit()

# calculate additional needed variables
preprocessor = PreProcessor()
preprocessor.process(events,
  leptongenvariables=['isPrompt'],
  leptonvariables=[
    'jetPtRatio',
    'jetBTagDeepFlavor'
  ],
  topmvavariable='mvaTOP', topmvaversion='ULv1'
)

# do object and event selection
selected_muons = events.Muon[muonselection(events.Muon, selectionid='run2ul_loose')]
selected_electrons = events.Electron[
  ( (electronselection(events.Electron, selectionid='run2ul_loose'))
    & (clean_electrons_from_muons(events.Electron, selected_muons)) )
]
selected_events = events[
  (ak.num(selected_muons.pt) + ak.num(selected_electrons.pt) >= 2)
]

# write to a file
writer = NanoEventsWriter()
writer.write( selected_events, args.outputfile, drop='default' )

# print size of input and output file
print('Size of input file:')
os.system('ls -lh {}'.format(args.inputfile))
print('Size of output file:')
os.system('ls -lh {}'.format(args.outputfile))
