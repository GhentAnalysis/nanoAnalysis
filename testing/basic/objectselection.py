#####################################################
# Basic test script for performing object selection #
#####################################################
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

# check number of electrons
nprintevents = 10
print('Number of electrons:')
print(ak.num(events.Electron, axis=1)[:nprintevents])

# do object selection
selected_electrons = events.Electron[
  (events.Electron.pt > 30) 
  & (abs(events.Electron.eta) < 2.1)]

# check type of created objects
print('Info on selected_electrons:')
print(selected_electrons)
print(type(selected_electrons))

# check number of electrons after selection
print('Number of electrons:')
print(ak.num(selected_electrons, axis=1)[:nprintevents])

# alternative approach: do not create a separate array
# with selected electrons, only define a selection mask
electron_selection = ((events.Electron.pt > 30) 
                     & (abs(events.Electron.eta) < 2.1 ))
print('Alternative approach')
print(type(electron_selection))
print(ak.num(events.Electron[electron_selection], axis=1)[:nprintevents])
