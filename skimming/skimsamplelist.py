##################################
# Skim tuples from a sample list #
##################################

import sys
import os
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
import jobsubmission.condortools as ct
from jobsubmission.jobsettings import CMSSW_VERSION
import tools.argparsetools as apt
from samples.samplelisttools import readsamplelist
from tools.dastools import get_sample_files
from tools.listtools import makechunks


# parse arguments
parser = argparse.ArgumentParser('Skim samples from a list')
parser.add_argument('-i', '--samplelist', required=True, type=os.path.abspath)
parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
parser.add_argument('-n', '--nentries', type=int, default=-1)
parser.add_argument('-l', '--leptonselection', default=None)
parser.add_argument('-s', '--skimselection', default=None)
parser.add_argument('-d', '--dropbranches', default=None)
parser.add_argument('--selectfirst', default=False, action='store_true')
parser.add_argument('--twosteps', default=False, action='store_true')
parser.add_argument('--compressionlevel', default=9, type=int)
parser.add_argument('--files_per_job', default=10, type=int)
parser.add_argument('--walltime_hours', default=24, type=int)
parser.add_argument('--filemode', default='das', choices=['das','local'])
parser.add_argument('--inputdir', default=None, type=apt.path_or_none)
parser.add_argument('--proxy', default=None, type=apt.path_or_none)
parser.add_argument('--max_files_per_sample', default=-1, type=int)
parser.add_argument('--readmode', default='remote', choices=['remote','copy'])
parser.add_argument('--runmode', default='condor', choices=['condor','local'])
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# check if input directory exists
if args.inputdir is not None:
    if not os.path.exists(args.inputdir):
        raise Exception('ERROR: input directory {} does not exist.'.format(args.inputdir))

# check if sample list exist
if not os.path.exists(args.samplelist):
    raise Exception('ERROR: sample list {} does not exist.'.format(args.samplelist))

# check if output directory is empty and ask permission to clean it
if os.path.exists(args.outputdir):
    if not len(os.listdir(args.outputdir))==0:
        print('WARNING: output directory {} not empty.'.format(args.outputdir)
              +' Permission to clean it? (y/n)')
        go = input()
        if not go=='y': sys.exit()
        os.system('rm -r {}'.format(os.path.join(args.outputdir,'*')))
else: os.makedirs(args.outputdir)

# parse walltime
walltime = '{}:00:00'.format(args.walltime_hours)

# make a list of sample names
print('Extracting sample names from the provided sample list...')
sample_collection = readsamplelist(args.samplelist, doyear=False)
sample_names = [s.name for s in sample_collection.get_samples()]

# print sample names
print('Extracted {} sample names from the sample list:'.format(len(sample_names)))
for s in sample_names: print('  - {}'.format(s))

# get the files for each sample
print('Finding number of files to process...')
sample_files = {}
nfiles = []
for s in sample_names:
    maxfiles = None
    if args.max_files_per_sample > 0: maxfiles = args.max_files_per_sample
    this_sample_files = get_sample_files(s,
                      filemode=args.filemode,
                      maxfiles=maxfiles)
    # check if sample was found correctly
    allgood = True
    if len(this_sample_files)==0: allgood = False
    for this_sample_file in this_sample_files:
        if not this_sample_file.endswith('.root'): allgood = False
    if not allgood:
        msg = 'ERROR: something seems wrong with sample {}'.format(s)
        msg += ' query to retrieve files for this sample'
        msg += ' returned the following: {}'.format(this_sample_files)
        raise Exception(msg)
    # do printouts for checking
    print('sample: {} -> {} files found'.format(s,len(this_sample_files)))
    # add files to list
    sample_files[s] = this_sample_files
    nfiles.append(len(this_sample_files))
# estimate number of jobs and do printouts for checking
nfiles = sum(nfiles)
njobs = max(1,int(nfiles/args.files_per_job))
print('Found a total of {} files, which will result in approximately {} jobs.'.format(
    nfiles,njobs))
print('Continue with the submission? (y/n)')
go = input()
if not go=='y': sys.exit()

# make output directory for each sample
print('Making output directories...')
sample_output_directories = []
for sample_name in sample_names:
    _, shortname, version, _ = sample_name.split('/')
    output_directory = os.path.join( args.outputdir, 
      'ntuples_skimmed_{}_version_{}'.format( shortname, version ) )
    if not os.path.exists( output_directory ): os.makedirs( output_directory )
    sample_output_directories.append( output_directory )

# loop over samples and submit skimming jobs
print('Starting submission...')
cwd = os.getcwd()
itlist = zip(sample_names, sample_output_directories)
for sample_name, sample_output_directory in itlist:
    print('Now processing the following sample:')
    print('  {}'.format(sample_name))
    print('  Number of root files: {}'.format(len(sample_files[sample_name])))
    # split files in lists of files_per_job
    chunks = list(makechunks( sample_files[sample_name], args.files_per_job ))
    for chunk in chunks:
        # make the commands to execute for this chunk
        commands = []
        commands.append( 'cd {}'.format(cwd) )
        # loop over files in this chunk
        for f in chunk:
            # define output file
            output_file = f.split('/')[-1]
            output_file = os.path.join(sample_output_directory,output_file)
            # define command to skim this file
            # (leave intput file blank as it will be added later depending on readmode)
            skimcommand = 'python3 skimfile.py -o {}'.format(output_file)
            if args.nentries > 0: skimcommand += ' -n {}'.format(args.nentries)
            if args.leptonselection is not None: skimcommand += ' -l {}'.format(args.leptonselection)
            if args.skimselection is not None: skimcommand += ' -s {}'.format(args.skimselection)
            if args.dropbranches is not None: skimcommand += ' -d {}'.format(args.dropbranches)
            if args.selectfirst: skimcommand += ' --selectfirst'
            if args.twosteps: skimcommand += ' --twosteps'
            if args.compressionlevel is not None: skimcommand += ' --compressionlevel {}'.format(args.compressionlevel)
            thiscommands = []
            if args.readmode=='remote':
                # read remote file directly
                skimcommand += ' -i {}'.format(f)
                thiscommands.append(skimcommand)
            elif args.readmode=='copy':
                # copy remote file to local before running skimmer
                output_file_unskimmed = output_file.replace('.root','_raw.root')
                thiscommands.append('xrdcp {} {}'.format(f, output_file_unskimmed))
                skimcommand += ' -i {}'.format(output_file_unskimmed)
                thiscommands.append(skimcommand)
                thiscommands.append('rm -f {}'.format(output_file_unskimmed))
            for c in thiscommands: commands.append(c)
        # run in local
        if( args.runmode=='local' ):
            for cmd in commands: os.system(cmd)
        # submission via condor
        if( args.runmode=='condor' ): 
            ct.submitCommandsAsCondorJob( 
              'cjob_skimsamplelist', commands,
              cmssw_version=CMSSW_VERSION, proxy=args.proxy )
