###################################################################
# Runs mergehists.py for a number of predefined regions and years #
###################################################################

import sys
import os

topdir = sys.argv[1]

regions = ['auto']

years = ['2018']

npmodes = []
#npmodes.append( 'npfromsim' )
npmodes.append( 'npfromdata' )

cfmodes = []
#cfmodes.append( 'cffromsim' )
cfmodes.append( 'cffromdata' )

rename = 'processes/rename_processes.json'
renamemode = 'fast'

decorrelate = 'correlations/correlations.json'
decorrelatemode = 'fast'

selectmode = 'noselect'

doclip = True

runmode = 'condor'

cmds = []
for year in years:
  yeardir = os.path.join(topdir, year+'_binned')
  if 'auto' in regions: regions = os.listdir(yeardir)
  for region in regions:
    for npmode in npmodes:
      for cfmode in cfmodes:
        inputdir = os.path.join(yeardir, region)
        outputfile = os.path.join(yeardir, region, 'merged_{}_{}'.format(npmode,cfmode), 'merged.root')
        cmd = 'python mergehists.py'
        cmd += ' -d '+inputdir
        cmd += ' -o '+outputfile
        cmd += ' --npmode '+npmode
        cmd += ' --cfmode '+cfmode
        cmd += ' --split'
        if rename is not None:
          cmd += ' --rename '+rename
          cmd += ' --renamemode '+renamemode
        if decorrelate is not None:
          cmd += ' --decorrelate '+decorrelate
          cmd += ' --decorrelatemode '+decorrelatemode
          cmd += ' --decorrelateyear '+year
        cmd += ' --selectmode '+selectmode
        if doclip: cmd += ' --doclip'
        cmd += ' --runmode '+runmode
        cmds.append(cmd)

for cmd in cmds:
  print('executing '+cmd)
  os.system(cmd)
