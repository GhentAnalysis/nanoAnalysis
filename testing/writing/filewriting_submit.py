################################
# Submitter for filewriting.py #
################################

import sys
import os
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION

# input arguments:
# input arguments:
parser = argparse.ArgumentParser(description='Test file writing and re-reading')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('-r', '--runmode', choices=['condor','local'], default='condor')
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make the command
cmd = 'python3 filewriting.py'
cmd += ' -i {}'.format(args.inputfile)
cmd += ' -o {}'.format(args.outputfile)
if args.nentries>0: cmd += ' -n {}'.format(args.nentries)

# submit job
if args.runmode=='condor':
    ct.submitCommandAsCondorJob( 'cjob_filewriting', cmd,
      cmssw_version=CMSSW_VERSION )
else:
    os.system(cmd)
