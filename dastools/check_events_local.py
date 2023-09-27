################################################################
# Check if specified events are in a dataset on stored locally #
################################################################
# Note: basically the same as check_events_das but for locally accessible 

# imports
import sys
import os
import argparse
import numpy as np
import uproot


if __name__=='__main__':
  
    parser = argparse.ArgumentParser(description='Check if event is in dataset')
    parser.add_argument('-d', '--datadir', required=True,
                        help='Directory holding one or more nanoAOD files')
    parser.add_argument('-e', '--events', required=True)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # read queried events
    events = []
    print('Reading events to check')
    with open(args.events, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip(' \t\n')
        if ':' in line: lineparts = line.split(':')
        else: lineparts = line.split()
        if not len(lineparts)==3:
            raise Exception('ERROR: events file seems to have unexpected format.')
        events.append((int(lineparts[0]), int(lineparts[1]), int(lineparts[2])))
    print('Found {} events to check'.format(len(events)))

    # find files
    rfiles = [f for f in os.listdir(args.datadir) if f.endswith('.root')]
    rfiles = [os.path.join(args.datadir,f) for f in rfiles]
    print('Found {} files'.format(len(rfiles)))

    # read all files
    runs = []
    lumis = []
    eventnbs = []
    print('Reading all files')
    for i, rfile in enumerate(rfiles):
        print('  reading file {} out of {}'.format(i+1, len(rfiles)))
        with uproot.open(rfile) as f:
            tree = f['Events']
            runs.append( tree['run'].array(library='np') )
            lumis.append( tree['luminosityBlock'].array(library='np') )
            eventnbs.append( tree['event'].array(library='np') )
    runs = np.concatenate(runs)
    lumis = np.concatenate(lumis)
    eventnbs = np.concatenate(eventnbs)

    # loop over queried events
    foundevents = []
    notfoundevents = []
    for event in events:
        print('Event: {}'.format(event))
        found = np.any((runs==event[0]) & (lumis==event[1]) & (eventnbs==event[2]))
        if found:
            print('-> found!')
            foundevents.append(event)
        else:
            print('-> not found!')
            notfoundevents.append(event)

    # print summary
    print('Summary:')
    print('Found {} out of {} events in {}'.format(len(foundevents),len(events),args.datadir))
