import sys
import os
# import framework modules
sys.path.append(os.path.abspath('../../'))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION


if __name__=='__main__':

  runmode = 'condor'

  files = [
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaod/TTWJetsToLNu-RunIISummer20UL18-nanoAOD-filtered-skimmed.root',
    # 'nevents': 1e5,
    # 'tag': 'ttw'},
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/Data_Run2018A.root',
    # 'nevents': 1e6,
    # 'tag': 'data'},
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/DoubleMuon_Run2018D.root',
    # 'nevents': -1,
    # 'tag': 'doublemuon'},
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/SingleMuon_Run2018D.root',
    # 'nevents': -1,
    # 'tag': 'singlemuon'},
    {'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/SingleMuon_Run2018C.root',
      'nevents': -1,
      'tag': 'singlemuon'},
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/EGamma_Run2018D.root',
    # 'nevents': -1,
    # 'tag': 'egamma'},
    #{'inputfile': '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged/MuonEG_Run2018D.root',
    # 'nevents': -1,
    # 'tag': 'muoneg'}
  ]

  #selection_types = ['tight', 'fakerate', 'chargeflips']
  selection_types = ['fakerate']

  for f in files:
    for t in selection_types:
      cmds = []
      outputname = 'output_cutflow_{}_{}'.format(f['tag'], t)
      # fill the cutflow histogram
      cmd = 'python3 cutflow_fill.py'
      cmd += ' -i {}'.format(f['inputfile'])
      cmd += ' -s {}'.format('signalregion_dilepton_inclusive')
      cmd += ' -t {}'.format(t)
      cmd += ' -o {}'.format('{}.root'.format(outputname))
      cmd += ' -n {}'.format(int(f['nevents']))
      cmd += ' --skimmed'
      cmds.append(cmd)
      # make a plot
      cmd = 'python3 cutflow_plot.py'
      cmd += ' -i {}'.format('{}.root'.format(outputname))
      cmd += ' -o {}'.format('{}_plots'.format(outputname))
      #cmds.append(cmd)
      # run or submit the commands
      if runmode=='local':
        for cmd in cmds: os.system(cmd)
      elif runmode=='condor':
        ct.submitCommandsAsCondorJob('cjob_cutflow', cmds, cmssw_version=CMSSW_VERSION)
