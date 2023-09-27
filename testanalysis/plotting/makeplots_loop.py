##################################################################
# Runs makeplots.py for a number of predefined regions and years #
##################################################################

# import python modules
import sys
import os
# import framework modules
sys.path.append(os.path.abspath('../../jobsubmission'))
import condortools as ct
from jobsettings import CMSSW_VERSION

inputdir = sys.argv[1]
runmode = 'local'

regions = []
for r in ['signalregion_dilepton_inclusive']: regions.append(r)
#for r in ['signalregion_trilepton']: regions.append(r)
#for r in ['trileptoncontrolregion']: regions.append(r)
#for r in ['fourleptoncontrolregion']: regions.append(r)
for r in ['npcontrolregion_met_dilepton_inclusive']: regions.append(r)
#for r in ['npcontrolregion_lownjets_dilepton_inclusive']: regions.append(r)
#for r in ['cfcontrolregion_inclusivejets']: regions.append(r)
#for r in ['cfcontrolregion_highnjets']: regions.append(r)

years = []
#years = ['2016PreVFP','2016PostVFP','2017','2018']
years = ['2018']

npmodes = ['npfromsim','npfromdata']
cfmodes = ['cffromsim','cffromdata']

dolog = True

unblind = True

variables = '../variables/variables_test.json'

colormap = 'ttw'

signals = 'TTW'

cmds = []
for year in years:
  for npmode in npmodes:
    for cfmode in cfmodes:
      for region in regions:
        subdir = os.path.join(year+'_binned', 'merged_{}_{}'.format(npmode,cfmode))
        inputfile = os.path.join(inputdir, subdir, 'merged.root')
        if not os.path.exists(inputfile):
          print('WARNING: input file {} does not exist; continuing...'.format(inputfile))
          continue
        thisoutputdir = '{}_{}_{}_{}'.format(year,region,npmode,cfmode)
        thisoutputdir = os.path.join(inputdir, subdir, 'plots', thisoutputdir)
        cmd = 'python makeplots.py'
        cmd += ' --inputfile '+inputfile
        cmd += ' --year '+year
        cmd += ' --region '+region
        cmd += ' --variables '+variables
        cmd += ' --outputdir '+thisoutputdir
        if unblind: cmd += ' --unblind'
        if dolog: cmd += ' --dolog'
        cmd += ' --colormap '+colormap
        cmd += ' --signals '+signals
        if runmode=='local':
          print('executing '+cmd)
          os.system(cmd)
        elif runmode=='condor':
          print('submitting '+cmd)
          cmds.append(cmd)
        else: raise Exception('ERROR: runmode "{}" not recognized'.format(runmode))

if runmode=='condor':
  ct.submitCommandsAsCondorCluster('cjob_makeplots', cmds,
                                    cmssw_version=CMSSW_VERSION)
