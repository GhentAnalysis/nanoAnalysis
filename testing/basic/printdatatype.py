########################################################################
# Basic test script for reading a nanoAOD file into a NanoEvents array #
########################################################################
# See more info in the official coffea documentation:
# https://coffeateam.github.io/coffea/notebooks/nanoevents.html

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
    entry_stop=1,
    schemaclass=NanoAODSchema,
).events()
print('Basic properties of "events" array:')
print(events)

# print out some data types
print(events.Electron.pt)
print(type(events.Electron.pt))
print(events.Electron.pt.type)
