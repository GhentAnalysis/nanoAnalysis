########################################
# Run binner.py for a batch of samples #
########################################

# import python modules
import sys
import os
import argparse
from pathlib import Path
# import framework modules
sys.path.append(str(Path(__file__).parents[1]))
from samples.samplelisttools import readsamplelist
import tools.argparsetools as apt
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION


if __name__=='__main__':

  # parse arguments
  parser = argparse.ArgumentParser(description='Perform event binning')
  parser.add_argument('-i', '--inputdir', required=True, type=os.path.abspath)
  parser.add_argument('-l', '--samplelist', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
  parser.add_argument('-v', '--variables', required=True, type=os.path.abspath)
  parser.add_argument('-s', '--eventselection', default=['auto'], nargs='+')
  parser.add_argument('-t', '--selectiontype', default=['auto'], nargs='+')
  parser.add_argument('-n', '--nevents', default=0, type=int)
  parser.add_argument('--systematics', default=['auto'], nargs='+')
  parser.add_argument('--lumi', default=1, type=float)
  parser.add_argument('--split', default=False, action='store_true')
  parser.add_argument('--runmode', default='condor', choices=['condor','local'])
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

  # argument checks and parsing
  if not os.path.exists(args.inputdir):
    raise Exception('ERROR: input directory {} does not exist.'.format(args.inputdir))
  if not os.path.exists(args.samplelist):
    raise Exception('ERROR: sample list {} does not exist.'.format(args.samplelist))
  if not os.path.exists(args.variables):
    raise Exception('ERROR: variable file {} does not exist.'.format(args.variables))

  # make output directory
  if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

  # check samples
  samples = readsamplelist( args.samplelist, sampledir=args.inputdir, doyear=True )
  year = samples.samples[0].year
  nsamples = samples.number()
  print('Found {} samples:'.format(nsamples))
  print(samples)

  # loop over input files and submit jobs
  cmds = []
  for i, sample in enumerate(samples.samples):
    # define input and output file
    inputfile = sample.path
    outputfile = os.path.basename(inputfile).replace('.root', '_binned.root')
    outputfile = os.path.join(args.outputdir, outputfile)
    # find cross-section
    xsec = sample.xsec
    # find process name
    process = sample.process
    # make the command
    cmd = 'python3 binner.py'
    cmd += ' -i {}'.format(inputfile)
    cmd += ' -o {}'.format(outputfile)
    cmd += ' -v {}'.format(args.variables)
    cmd += ' -s {}'.format(' '.join(args.eventselection))
    cmd += ' -t {}'.format(' '.join(args.selectiontype))
    if len(args.systematics) > 0: cmd += ' --systematics {}'.format(' '.join(args.systematics))
    cmd += ' --xsec {}'.format(xsec)
    cmd += ' --lumi {}'.format(args.lumi)
    cmd += ' --process {}'.format(process)
    if args.split: cmd += ' --split'
    if args.nevents > 0: cmd += ' -n {}'.format(args.nevents)
    cmds.append(cmd)

  # submit the jobs
  if args.runmode=='local':
    for cmd in cmds: os.system(cmd)
  elif args.runmode=='condor':
    ct.submitCommandsAsCondorCluster( 'cjob_binner', cmds,
                                      cmssw_version=CMSSW_VERSION )
