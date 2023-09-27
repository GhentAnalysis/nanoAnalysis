##################################################################
# Runs eventloop.py for a number of predefined regions and years #
##################################################################

import os
import sys

# command line arguments

outputdir = sys.argv[1]

# hard coded arguments

regions = []
for r in ['signalregion_dilepton_inclusive']: regions.append(r)
for r in ['signalregion_trilepton']: regions.append(r)
for r in ['trileptoncontrolregion']: regions.append(r)
for r in ['fourleptoncontrolregion']: regions.append(r)
for r in ['npcontrolregion_met_dilepton_inclusive']: regions.append(r)
for r in ['npcontrolregion_lownjets_dilepton_inclusive']: regions.append(r)
for r in ['cfcontrolregion_inclusivejets']: regions.append(r)
for r in ['cfcontrolregion_highnjets']: regions.append(r)

years = ['2018']

selectiontypes = []
selectiontypes.append('tight')
selectiontypes.append('fakerate')
selectiontypes.append('chargeflips')
selectiontypes.append('irreducible')

dtypes = []
dtypes.append('sim')
dtypes.append('data')

frdir = '../data/fakerates/fakeRateMaps_v20220912_tttt'
cfdir = '../data/chargefliprates/chargeFlipMaps_v20221109'

samplelistdir = 'samplelists'
samplelistbase = {'sim': 'samplelist_ttw_{}_sim.txt',
                  'data': 'samplelist_ttw_{}_datasplit.txt'}

bdtfile = None

nevents = 1e5

skimmed = True

systematics = ['all']

runmode = 'condor'

for year in years:
  for dtype in dtypes:
    # set correct input directory
    inputdir = '/pnfs/iihe/cms/store/user/llambrec/nanoaodskims_merged'
    # set correct output directory
    thisoutputdir = os.path.join(outputdir, '{}_{}_trees'.format(year, dtype))
    # set correct sample list
    samplelist = os.path.join(samplelistdir,samplelistbase[dtype].format(year))
    # make basic command
    cmd = 'python3 eventloop_batch.py'
    cmd += ' -i ' + inputdir
    cmd += ' -l ' + samplelist
    cmd += ' -o ' + thisoutputdir
    cmd += ' -s {}'.format(' '.join(regions))
    cmd += ' -t {}'.format(' '.join(selectiontypes))
    cmd += ' --systematics {}'.format(' '.join(systematics))
    cmd += ' --frdir ' + frdir
    cmd += ' --cfdir ' + cfdir
    if  skimmed: cmd += ' --skimmed'
    cmd += ' --runmode ' + runmode
    if( nevents is not None and nevents > 0 ): cmd += ' -n ' + str(int(nevents))
    if bdtfile is not None: cmd += ' --bdt ' + bdtfile
    print('executing '+cmd)
    os.system(cmd)
