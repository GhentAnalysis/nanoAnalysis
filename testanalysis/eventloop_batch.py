###########################################
# Run eventloop.py for a batch of samples #
###########################################

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
# import local modules
sys.path.append(os.path.abspath('systematics'))
from systematics_type import systematics_type


if __name__=='__main__':

  # parse arguments
  parser = argparse.ArgumentParser(description='Perform TTW event selection and calculate analysis variables')
  parser.add_argument('-i', '--inputdir', required=True, type=os.path.abspath)
  parser.add_argument('-l', '--samplelist', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
  parser.add_argument('-s', '--eventselection', required=True, nargs='+')
  parser.add_argument('-t', '--selectiontype', default=['tight'], nargs='+')
  parser.add_argument('-n', '--nevents', default=0, type=int)
  parser.add_argument('--systematics', default=[], nargs='+')
  parser.add_argument('--frdir', default=None, type=apt.path_or_none)
  parser.add_argument('--cfdir', default=None, type=apt.path_or_none)
  parser.add_argument('--bdt', default=None, type=apt.path_or_none)
  parser.add_argument('--skimmed', default=False, action='store_true')
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
  if os.path.exists(args.outputdir):
    print('WARNING: output directory already exists. Clean it? (y/n)')
    go=input()
    if not go=='y': sys.exit()
    os.system('rm -r '+args.outputdir)

  # check systematics
  if( len(args.systematics)==1 and args.systematics[0]=='all' ): pass
  else:
    for systematic in args.systematics:
      if systematic not in systematics_type.keys():
        raise Exception('ERROR: systematic {} not recognized.'.format(systematic))

  # make output directory
  os.makedirs(args.outputdir)

  # check samples
  samples = readsamplelist( args.samplelist, sampledir=args.inputdir, doyear=True )
  year = samples.samples[0].year
  nsamples = samples.number()
  print('Found {} samples:'.format(nsamples))
  print(samples)

  # set and check fake rate maps
  muonfrmap = None
  electronfrmap = None
  if 'fakerate' in args.selectiontype:
    if args.frdir is None:
      raise Exception('ERROR: fake rate dir must be specified for selection type fakerate.')
    muonfrmap = os.path.join(args.frdir,'fakeRateMap_data_muon_'+year+'_mT.root')
    electronfrmap = os.path.join(args.frdir,'fakeRateMap_data_electron_'+year+'_mT.root')
    if not os.path.exists(muonfrmap):
      raise Exception('ERROR: fake rate map {} does not exist'.format(muonfrmap))
    if not os.path.exists(electronfrmap):
      raise Exception('ERROR: fake rate map {} does not exist'.format(electronfrmap))

  # set and check charge flip maps
  electroncfmap = None
  if 'chargeflips' in args.selectiontype:
    if args.cfdir is None:
      raise Exception('ERROR: charge flip dir must be specified for selection type chargeflips.')
    electroncfmap = os.path.join(args.cfdir,'chargeFlipMap_MC_electron_'+year+'.root')
    if not os.path.exists(electronfrmap):
      raise Exception('ERROR: fake rate map {} does not exist'.format(electronfrmap))

  # check bdt weight file
  if( args.bdt is not None ):
    if not os.path.exists(args.bdt):
      raise Exception('ERROR: BDT file {} does not exist'.format(args.bdt))

  # loop over input files and submit jobs
  cmds = []
  for i, sample in enumerate(samples.samples):
    # define input and output file
    inputfile = sample.path
    outputfile = os.path.join(args.outputdir, os.path.basename(inputfile))
    # make the command
    cmd = 'python3 eventloop.py'
    cmd += ' -i {}'.format(inputfile)
    cmd += ' -o {}'.format(outputfile)
    cmd += ' -s {}'.format(' '.join(args.eventselection))
    cmd += ' -t {}'.format(' '.join(args.selectiontype))
    if args.nevents > 0: cmd += ' -n {}'.format(args.nevents)
    if len(args.systematics) > 0: cmd += ' --systematics {}'.format(' '.join(args.systematics))
    if muonfrmap is not None: cmd += ' --mufrmap {}'.format(muonfrmap)
    if electronfrmap is not None: cmd += ' --elfrmap {}'.format(electronfrmap)
    if electroncfmap is not None: cmd += ' --elcfmap {}'.format(electroncfmap)
    if args.bdt is not None: cmd += ' --bdt {}'.format(args.bdt)
    if args.skimmed: cmd += ' --skimmed'
    cmds.append(cmd)

  # submit the jobs
  if args.runmode=='local':
    for cmd in cmds: os.system(cmd)
  elif args.runmode=='condor':
    ct.submitCommandsAsCondorCluster( 'cjob_eventloop', cmds,
                                      cmssw_version=CMSSW_VERSION )
