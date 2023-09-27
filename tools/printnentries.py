####################################################
# Print number of entries in given NanoAOD file(s) #
####################################################

# imports
import sys
import os
import awkward as ak
import uproot

# input arguments: arbitrary number of files
inputfiles = sys.argv[1:]

# loop over input files
for inputfile in inputfiles:

  # open the file with uproot
  with uproot.open(inputfile) as f:
  
    # get the number of entries in the Events tree
    events = f['Events']
    nentries = events.num_entries

    # do printouts
    print('File: {}: {} entries'.format(inputfile, nentries))
