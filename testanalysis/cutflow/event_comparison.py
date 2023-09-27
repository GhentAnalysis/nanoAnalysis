#########################################################
# Temporary script for comparing lists of event numbers #
#########################################################

import sys
import os
import numpy as np
import uproot

if __name__=='__main__':

  # command line args
  file1 = sys.argv[1]
  file2 = sys.argv[2]
  file3 = sys.argv[3]
  rootfile = sys.argv[4]

  # read lists of event numbers
  files = [file1, file2, file3]
  eventids = []
  for f in files:
    print('Reading {}'.format(f))
    with open(f, 'r') as s:
      lines = s.readlines()
    thiseventids = []
    for line in lines:
        line = line.strip(' \t\n')
        lineparts = line.split()
        #eventid = int(int(lineparts[0]) + 1e6*int(lineparts[1]) + 1e10*int(lineparts[2])) # use full number
        #eventid = int(lineparts[2]) # ignore run and LS number
        eventid = line # use str representation
        thiseventids.append(eventid)
    thiseventids = np.array(thiseventids)
    eventids.append(thiseventids)
    print('Found {} events'.format(len(thiseventids)))
  with uproot.open(rootfile) as f:
    print('Reading {}'.format(rootfile))
    tree = f['Events']
    runs = tree['run'].array(library='np')
    lumis = tree['luminosityBlock'].array(library='np')
    events = tree['event'].array(library='np')
  print('Found {} events'.format(len(events)))
  print('Formatting')
  # use full number
  #thiseventids = np.vstack((runs, lumis*1e6, events*1e10))
  #thiseventids = np.sum(thiseventids, axis=0)
  # ignore run and LS number
  #thiseventids = events
  # use str representation
  runs = np.char.add(runs.astype(str), ' ')
  lumis = np.char.add(lumis.astype(str), ' ')
  thiseventids = np.char.add(runs, lumis)
  thiseventids = np.char.add(thiseventids, events.astype(str))
  eventids.append(thiseventids)

  # find unique events
  print('Calculating differences')
  unique_events_1 = np.setdiff1d(eventids[0], eventids[1])
  unique_events_2 = np.setdiff1d(unique_events_1, eventids[2])
  unique_events_3 = np.setdiff1d(unique_events_2, eventids[3])

  # printouts
  print('Events in 1: {}'.format(len(eventids[0])))
  print('Events in 2: {}'.format(len(eventids[1])))
  print('Events in 1 and 2: {}'.format(len(np.intersect1d(eventids[0],eventids[1]))))
  print('Events in 1 but not in 2: {}'.format(len(unique_events_1)))
  print('Events in 2 but not in 1: {}'.format(len(np.setdiff1d(eventids[1], eventids[0]))))
  print('Events in 1 but not in 2 or 3: {}'.format(len(unique_events_2)))
  print('Events in 1 but not in 2 or 3 or 4: {}'.format(len(unique_events_3)))

  # write missing events
  with open('temp_missing.txt', 'w') as f:
    for eventid in unique_events_3:
      f.write('{}\n'.format(eventid))
