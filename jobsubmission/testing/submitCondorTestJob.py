################################################
# A tester for the condor job submission tools #
################################################
# How to use?
#   Simply run "python submitCondorTestJob.py" without arguments.
#   To first test in local, run "python submitCondorTestJob.py local".

import sys
import os
import time
import datetime
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION


def some_task():
    ### perform a basic task that takes some time

    print('###starting###')
    for i in range(5):
        ctime = datetime.datetime.now().strftime('%H:%M:%S')
        print('Current time: {}'.format(ctime))
        time.sleep(5)
    print('###done###')


if __name__=='__main__':

    nargs = len(sys.argv)
    if nargs==1:
        cmd = 'python submitCondorTestJob.py local'
        ct.submitCommandAsCondorJob('cjob_submitCondorTestJob', cmd,
            cmssw_version=CMSSW_VERSION)
    elif nargs==2:
        if sys.argv[1]=='local':
            some_task()
        else:
            raise Exception('ERROR: argument not recognized:'
                            +' {}'.format(sys.argv[1]))
    else:
        raise Exception('ERROR: unrecognized number of arguments:'
                            +' {} ({})'.format(sys.argv,nargs))
