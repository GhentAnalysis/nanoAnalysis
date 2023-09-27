###############################################################
# Runs binner.py for a number of predefined regions and years #
###############################################################

import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
from constants.luminosities import lumidict

inputdir = sys.argv[1]

regions = ['auto']
selectiontypes = ['auto']
systematics = ['auto']

years = ['2018']

dtypes = []
dtypes.append('sim')
dtypes.append('data')

samplelistdir = 'samplelists'
samplelistbase = {'sim': 'samplelist_ttw_{}_sim.txt',
                  'data': 'samplelist_ttw_{}_datamerged.txt'}

variables = 'variables/variables_test.json'

nevents = None

runmode = 'condor'

for year in years:
  for dtype in dtypes:
    # set correct input directory
    thisinputdir = os.path.join(inputdir, '{}_{}_trees'.format(year, dtype))
    # set output directory
    thisoutputdir = thisinputdir.replace('_trees', '_binned')
    # set correct sample list
    samplelist = os.path.join(samplelistdir,samplelistbase[dtype].format(year))
    # make basic command
    cmd = 'python3 binner_batch.py'
    cmd += ' -i ' + thisinputdir
    cmd += ' -l ' + samplelist
    cmd += ' -o ' + thisoutputdir
    cmd += ' -v ' + variables
    cmd += ' -s {}'.format(' '.join(regions))
    cmd += ' -t {}'.format(' '.join(selectiontypes))
    cmd += ' --lumi {}'.format(lumidict[year])
    cmd += ' --split'
    cmd += ' --runmode ' + runmode
    if( nevents is not None and nevents > 0 ): cmd += ' -n ' + str(int(nevents))
    print('executing '+cmd)
    os.system(cmd)
