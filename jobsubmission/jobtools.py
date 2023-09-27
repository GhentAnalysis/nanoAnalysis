########################################
# Other tools for handling condor jobs #
########################################

import sys
import os
import subprocess


def get_running_jobs(allusers=False, jobtag=None):
    # get output of condor_q command
    cmd = 'condor_q'
    cmdargs = []
    if allusers: cmdargs.append('-all')
    cmdres = subprocess.run([cmd]+cmdargs, stdout=subprocess.PIPE)
    cmdout = cmdres.stdout.decode('utf-8')
    lines = cmdout.split('\n')
    # find header
    headeridx = -1
    for i,line in enumerate(lines):
        if line.startswith('OWNER'):
            headeridx = i
            break
    if headeridx < 0:
        raise Exception('ERROR: could not find header index in condor_q output.')
    # find footer
    footeridx = -1
    for i,line in enumerate(lines):
        if line.startswith('Total for query:'):
            footeridx = i
            break
    if( footeridx < 0 or footeridx <= headeridx ):
        raise Exception('ERROR: could not find footer index in condor_q output.')
    # find all jobs between header and footer
    joblines = lines[headeridx+1:footeridx]
    joblines = [l for l in joblines if len(l)>0]
    # select jobs if requested
    if jobtag is not None: joblines = [l for l in joblines if jobtag in l]
    # return the result
    return joblines

def jobs_are_running(verbose=False, **kwargs):
    if verbose: print('Now running condor_q to check running jobs...')
    joblines = get_running_jobs(**kwargs)
    if len(joblines)>0:
        if verbose:
            print('Found following jobs ({}):'.format(len(joblines)))
            for j in joblines: print(j)
        return True
    if verbose: print('No jobs found that meet search criteria.')
    return False
