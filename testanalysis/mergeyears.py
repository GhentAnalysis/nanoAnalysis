###############################################################################
# Merge the output of eventbinner over several years to make full run 2 plots #
###############################################################################
# Note: naive implementation of simply adding nominal histograms together;
#       taking into account systematics and correlations to be implemented later
#       using combine tools.
# Note: should be run after mergehists_submit.py.

import sys
import os
import argparse

if __name__=='__main__':

  # parse arguments
  parser = argparse.ArgumentParser(description='Merge years')
  parser.add_argument('--directory', required=True, type=os.path.abspath)
  parser.add_argument('--filemode', default='combined', choices=['combined','split'])
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))
  
  # fixed arguments (for now)
  years = ['2016PreVFP', '2016PostVFP', '2017', '2018']
  npmodes = ['npfromdata', 'npfromdatasplit'] #['npfromsim', 'npfromdata', 'npfromdatasplit']
  cfmodes = ['cffromdata'] #['cffromsim', 'cffromdata']

  # find regions (if applicable)
  regions = ['']
  if args.filemode=='split': 
    regions = os.listdir(os.path.join(args.directory,years[0]))

  # loop over combinations
  for region in regions:
    for npmode in npmodes:
      for cfmode in cfmodes:
        # check all required files
        allfiles = True
        hfiles = []
        for year in years:
          f = os.path.join(args.directory,year,region,'merged_{}_{}'.format(npmode,cfmode),'merged.root')
          if not os.path.exists(f):
            msg = 'ERROR: file {} not does not exist (for year {});'.format(f,year)
            msg += ' skipping merging for npmode {} cfmode {}.'.format(npmode,cfmode)
            print(msg)
            allfiles = False
          else: hfiles.append(f)
        if not allfiles: continue
        # define output file
        outdir = os.path.join(args.directory,'run2',region,'merged_{}_{}'.format(npmode,cfmode))
        if not os.path.exists(outdir): os.makedirs(outdir)
        outf = os.path.join(outdir,'merged.root')
        # do hadd
        cmd = 'hadd -f {}'.format(outf)
        for f in hfiles: cmd += ' {}'.format(f)
        os.system(cmd)
