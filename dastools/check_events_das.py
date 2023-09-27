#####################################################
# Check if specified events are in a dataset on DAS #
#####################################################
# Note:
#   based on this implementation:
#   https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/Utilities/scripts/edmPickEvents.py
# Note:
#   DAS only allows to retrieve the file for a given lumisection,
#   next one needs to open the file and check if the event number is in there;
#   this step is format/tier-specific (i.e. NanoAOD in this case)

import sys
import os
import argparse
import numpy as np
import uproot

def get_event_file(datasetname, event):
    dasquery = 'file dataset={} run={} lumi={}'.format(datasetname, event[0], event[1])
    dascmd = "dasgoclient -query '{}'".format(dasquery)
    dasstdout = os.popen(dascmd).read()
    if isinstance(dasstdout, list):
        raise Exception('Found list, need to catch this exception')
    f = dasstdout.strip(' \t\n')
    if len(f)==0: return None
    return f

 
if __name__=='__main__':
  
    parser = argparse.ArgumentParser(description='Check if event is in dataset')
    parser.add_argument('-d', '--datasetname', required=True,
                        help='Name of the dataset on DAS')
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

    # find all unique lumisections
    alllumis = []
    for event in events:
        lumi = (event[0], event[1])
        if lumi in alllumis: continue
        alllumis.append(lumi)
    print('Found {} lumisections'.format(len(alllumis)))

    # make a dictionary matching lumis to files
    print('Retrieving correct file for each lumisection')
    lumi_to_file_dict = {}
    for i, lumi in enumerate(alllumis):
        print('  processing lumi {} out of {}'.format(i+1, len(alllumis)))
        lumi_to_file_dict[lumi] = get_event_file(args.datasetname, lumi)
    allfiles = []
    for val in lumi_to_file_dict.values():
        if val is None: continue
        if val in allfiles: continue
        allfiles.append(val)
    print('Found {} files to read'.format(len(allfiles)))

    # make a dictionary matching queried events to files using both above
    eventdict = {}
    for event in events:
        lumi = (event[0], event[1])
        eventdict[event] = lumi_to_file_dict[lumi]

    # do printout of event to file matching
    for event, rootfile in eventdict.items():
        print('  {} -> {}'.format(event, rootfile))

    # make a dictionary matching files to available events
    filedict = {}
    print('Reading all files')
    for i, rootfile in enumerate(allfiles):
        print('  reading file {} out of {}'.format(i+1, len(allfiles)))
        remotefile = 'root://cms-xrd-global.cern.ch//{}'.format(rootfile)
        uproot.open.defaults["xrootd_handler"] = uproot.MultithreadedXRootDSource
        with uproot.open(remotefile) as f:
            tree = f['Events']
            runs = tree['run'].array(library='np')
            lumis = tree['luminosityBlock'].array(library='np')
            eventnbs = tree['event'].array(library='np')
        filedict[rootfile] = {'runs': runs, 'lumis': lumis, 'eventnbs': eventnbs}

    # loop over queried events
    foundevents = []
    notfoundevents = []
    for event in events:
        print('Event: {}'.format(event))
        rootfile = eventdict[event]
        if rootfile is None:
            print('-> not found!')
            notfoundevents.append(event)
            continue
        runs = filedict[rootfile]['runs']
        lumis = filedict[rootfile]['lumis']
        eventnbs = filedict[rootfile]['eventnbs']
        mask = ((runs==event[0]) & (lumis==event[1]))
        eventnbs = eventnbs[mask]
        if event[2] in eventnbs:
            print('-> found!')
            foundevents.append(event)
        else:
            print('-> not found!')
            notfoundevents.append(event)

    # print summary
    print('Summary:')
    print('Found {} out of {} events in {}'.format(len(foundevents),len(events),args.datasetname))
