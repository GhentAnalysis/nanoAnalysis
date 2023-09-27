################################################
# Run mergedatatrees over a full output folder #
################################################

import sys
import os
from pathlib import Path
# import framework modules
sys.path.append(str(Path(__file__).parents[1]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION


if __name__=='__main__':

  inputdir = sys.argv[1]
  runmode = 'condor'
  years = ['2018']
 
  # loop over years
  cmds = []
  for year in years:
    thisinputdir = os.path.join(inputdir, '{}_data_trees'.format(year))
    if not os.path.exists(thisinputdir):
      msg = 'WARNING: folder {} does not exist, skipping...'.format(thisinputdir)
      print(msg)
      continue
    # find all files
    dfiles = [f for f in os.listdir(thisinputdir) if f.endswith('.root')]
    # group files to merge
    eras = sorted(list(set([f.split('_')[-1].replace('.root','') for f in dfiles])))
    iofiles = {}
    for era in eras:
      ofile = 'Data_{}.root'.format(era)
      ifiles = [f for f in dfiles if era in f]
      iofiles[ofile] = ifiles
    # printouts for testing
    print('Will merge as follows:')
    for key, val in iofiles.items():
      print('  - {}: {}'.format(key, val))
    # loop over files to merge
    for ofile, ifiles in iofiles.items():
      # make command
      cmd = 'python3 mergedatatrees.py'
      cmd += ' -o {}'.format(os.path.join(thisinputdir,ofile))
      cmd += ' -i {}'.format(' '.join([os.path.join(thisinputdir,ifile) 
        for ifile in ifiles]))
      cmd += ' -f -v'
      cmds.append(cmd)

  # run the commands or submit the jobs
  if runmode=='local':
    for cmd in cmds: os.system(cmd)
  elif runmode=='condor':
    ct.submitCommandsAsCondorCluster( 'cjob_mergedatatrees', cmds,
                                      cmssw_version=CMSSW_VERSION )
