#######################################
# Submitter for makeroccurves_fill.py #
#######################################

import sys
import os
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION

# input arguments:
parser = argparse.ArgumentParser(description='Make ROC curves for lepton MVAs')
parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
parser.add_argument('-y', '--year', required=True)
# (to do: tools for extracting year automatically from file name)
parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('-r', '--runmode', choices=['condor', 'local'], default='condor')
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# make the command
cmd = 'python3 makeroccurves_fill.py'
cmd += ' -i {}'.format(args.inputfile)
cmd += ' -y {}'.format(args.year)
cmd += ' -o {}'.format(args.outputfile)
if args.nentries>0: cmd += ' -n {}'.format(args.nentries)

# submit job
if args.runmode=='condor':
    ct.submitCommandAsCondorJob( 'cjob_makeroccurves_fill', cmd,
      cmssw_version=CMSSW_VERSION )
else:
    os.system(cmd)
