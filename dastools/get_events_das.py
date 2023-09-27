#################################################################
# Get runs / lumisections / event numbers from a dataset on DAS #
#################################################################
# Note:
#   Works for runs and lumisections; but does not seem possible for event numbers.
#   (See e.g. here: https://github.com/dmwm/DAS/issues/4220.)
#   A workaround is to retrieve all files in a dataset (with DAS client),
#   and then opening all the files and reading their event numbers
#   (see also check_events_das.py), but this has not yet been implemented.

import sys
import os
import argparse

def get_runs_das(datasetname):
    # returns:
    # list of run numbers (in string format)
    dascmd = "dasgoclient -query 'run dataset={}' --limit 0".format(datasetname)
    dasstdout = os.popen(dascmd).read()
    runs = sorted([el.strip(' \t') for el in dasstdout.strip('\n').split('\n')])
    # format of runlumis is a list with strings of format '<run nb>'
    if len(runs)==1 and runs[0]=='': return []
    return runs

def get_lumis_das(datasetname):
    # returns:
    # dict of run numbers (in string format) to list of lumisections (in int format)
    dascmd = "dasgoclient -query 'run lumi dataset={}' --limit 0".format(datasetname)
    dasstdout = os.popen(dascmd).read()
    runlumis = sorted([el.strip(' \t') for el in dasstdout.strip('\n').split('\n')])
    # format of runlumis is a list with strings of format '<run nb> [<ls nbs>]'
    if len(runlumis)==1 and runlumis[0]=='': return []
    runsls = {}
    for runlumi in runlumis:
        run = str(int(runlumi.split(' ',1)[0]))
        lumis = runlumi.split(' ',1)[1]
        lumis = lumis.strip('[] ')
        lumis = lumis.split(',')
        lumis = [int(lumi) for lumi in lumis]
        lumis = set(lumis)
        # (note: the above is to remove duplicates that are sometimes observed;
        #  not sure where it is coming from or what it means...)
        lumis = sorted(list(lumis))
        runsls[run] = lumis
    return runsls


if __name__=='__main__':
  
    parser = argparse.ArgumentParser(description='Check available runs/lumis/events')
    parser.add_argument('-d', '--datasetname', required=True,
                        help='Name of the dataset on DAS')
    parser.add_argument('-l', '--level', default='run', choices=['run','lumi','event'])
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # get the requested info
    res = []
    if args.level=='run':
        res = get_runs_das(args.datasetname)
    elif args.level=='lumi':
        res = get_lumis_das(args.datasetname)
    elif args.level=='event':
        raise Exception('Not yet implemented')
    else:
        raise Exception('ERROR: level {} not recognized.'.format(args.level))

    # write to output file (to implement)
    print(res)
