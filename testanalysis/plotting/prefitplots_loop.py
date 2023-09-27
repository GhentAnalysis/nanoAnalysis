########################################################################################
# Simple submitter that runs makeplots.py for a number of predefined regions and years #
########################################################################################

import sys
import os
sys.path.append('../../jobsubmission')
import condortools as ct
#from jobSettings import CMSSW_VERSION
CMSSW_VERSION = '~/CMSSW_10_6_29'

inputdir = sys.argv[1]
runmode = 'condor'

regions = ['auto']

years = []
years.append('2018')

npmodes = ['npfromdata']
cfmodes = ['cffromdata']

unblind = True

dummysystematics = False

rawsystematics = False

dolog = True

variables = '../variables/variables_test.json' # single variables

colormap = 'ttw'

datatag = 'Data'

signals = ['TTW'] # for single variables

cmds = []
for year in years:
  yeardir = os.path.join(inputdir, year+'_binned')
  if 'auto' in regions: regions = os.listdir(yeardir)
  for region in regions:
    for npmode in npmodes:
      for cfmode in cfmodes:
          subdir = os.path.join(yeardir, region, 'merged_{}_{}'.format(npmode,cfmode))
          inputfile = os.path.join(inputdir, subdir, 'merged.root')
          if not os.path.exists(inputfile):
            print('WARNING: input file {} does not exist; continuing...'.format(inputfile))
            continue
          thisoutputdir = '{}_{}_{}_{}'.format(year,region,npmode,cfmode)
          if not unblind: thisoutputdir += '_blind'
          if rawsystematics: thisoutputdir += '_rawsystematics'
          if dummysystematics: thisoutputdir += '_dummysystematics'
          thisoutputdir = os.path.join(inputdir, subdir, 'plots', thisoutputdir)
          cmd = 'python prefitplots.py'
          cmd += ' --inputfile '+inputfile
          cmd += ' --year '+year
          cmd += ' --region '+region
          cmd += ' --processes all'
          cmd += ' --variables '+variables
          cmd += ' --outputdir '+thisoutputdir
          cmd += ' --datatag '+datatag
          cmd += ' --colormap '+colormap
          if unblind: cmd += ' --unblind'
          if rawsystematics: cmd += ' --rawsystematics'
          if dummysystematics: cmd += ' --dummysystematics'
          if dolog: cmd += ' --dolog'
          if signals is not None:
            cmd += ' --signals '+' '.join(signals)
          if runmode=='local':
            print('executing '+cmd)
            os.system(cmd)
          elif runmode=='condor':
            print('submitting '+cmd)
            cmds.append(cmd)
          else: raise Exception('ERROR: runmode "{}" not recognized'.format(runmode))

if runmode=='condor':
  ct.submitCommandsAsCondorCluster('cjob_prefitplots', cmds,
                                    cmssw_version=CMSSW_VERSION)
