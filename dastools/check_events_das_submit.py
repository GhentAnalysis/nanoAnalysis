#####################################
# Submitter for check_events_das.py #
#####################################
# Note:
#   Does not work yet...
#   Problem with DAS queries from inside a job...


# import python modules
import sys
import os
import argparse
from pathlib import Path
# import framework modules
sys.path.append(str(Path(__file__).parents[1]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser(description='Check if event is in dataset')
    parser.add_argument('-d', '--datasetname', required=True,
                        help='Name of the dataset on DAS')
    parser.add_argument('-e', '--events', required=True, type=os.path.abspath)
    parser.add_argument('-p', '--proxy', required=True, type=os.path.abspath)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # make the command
    cmd = 'python3 check_events_das.py'
    cmd += ' -d {}'.format(args.datasetname)
    cmd += ' -e {}'.format(args.events)
    print('Submitting following command:')
    print(cmd)

    # submit the command
    ct.submitCommandAsCondorJob( 'cjob_check_events_das', cmd,
                                 cmssw_version=CMSSW_VERSION,
                                 proxy=args.proxy )
