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
    schemaclass=NanoAODSchema,
).events()
print('Basic properties of "events" array:')
print(events)
print(type(events))

# print out all the fields
print('Contents of "events" array:')
for field in sorted(events.fields):
    print('  {}'.format(field))
    for subfield in sorted(events[field].fields):
        print('    {}'.format(subfield))
