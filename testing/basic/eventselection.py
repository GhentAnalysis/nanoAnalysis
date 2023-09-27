####################################################
# Basic test script for performing event selection #
####################################################
# Also including object selection, see objectselection.py
# for a more basic example.
# See more info in the cms open data workshop:
# https://cms-opendata-workshop.github.io/workshop2022-lesson-ttbarljetsanalysis/02-coffea-analysis/

# imports
import sys
import os
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# input arguments:
# - input file (expected to be in nanoAOD format)
inputfile = sys.argv[1]

# make NanoEvents array
events = NanoEventsFactory.from_root(
    inputfile,
    schemaclass=NanoAODSchema,
).events()

# define object selection
electron_selection = (
  (events.Electron.pt > 30) 
  & (abs(events.Electron.eta) < 2.1)
)
muon_selection = (
  (events.Muon.pt > 30)
  & (abs(events.Muon.eta) < 2.1)
)
jet_selection = (
  (events.Jet.pt > 30)
  & (abs(events.Jet.eta) < 2.4)
)

# define event selection
# (dummy selection with at least 1 electron, 1 muon and 3 jets)
event_selection = (
  (ak.num(events.Electron[electron_selection], axis=1) >= 1)
  & (ak.num(events.Muon[muon_selection], axis=1) >= 1)
  & (ak.num(events.Jet[jet_selection], axis=1) >= 3)
)

# check number of events before and after selection
print('Number of events before selection:')
print(ak.num(events.event, axis=0))
print('Number of events after selection:')
print(ak.num(events[event_selection].event, axis=0))

# check number of objects in selected events
print('Number of selected electrons, muons and jets in selected events:')
print(ak.num(events.Electron[electron_selection][event_selection], axis=1))
print(ak.num(events.Muon[muon_selection][event_selection], axis=1))
print(ak.num(events.Jet[jet_selection][event_selection], axis=1))
