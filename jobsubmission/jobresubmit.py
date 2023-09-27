#!/usr/bin/env python

#################################
# Resubmit stuck or failed jobs #
#################################
# Use case:
#   Use this tool for automated resubmission of failed jobs.
#   Schematically, this script does the following:
#   - Check the _err_ files for unfinished jobs (see jobCheck.py).
#   - Perform "condor_rm" on all jobs with errors.
#     (Might print errors if they were already killed before, but that is ok.)
#   - For all _err_ files with errors, find the corresponding executable
#     from the printouts in the corresponding _out_ file.
#     (The name of the executable is printed on the first line, see condorTools.py.)
#   - Modify the job description file for each of these
#     and call "condor_submit" on this file.
#   - Remove the original _err_, _out_ and _log_ files of the resubmitted jobs.
# Usage:
#   Run 'python jobResubmission.py -h' for a list of options.
#   You can run the script from this directory and specify the job directory in the args,
#   or alternatively you can run the script from the job directory
#   (using 'python [path]/jobResubmission.py) and leave the directory arg at its default.
#   If you have sourced the 'source.sh' script in the project's main directory,
#   you can simply run 'jobResubmission [+args]' from anywhere, 
#   without specifying "python" or the path to this script.
# Note: 
#   (So far) only for condor jobs, older qsub is not supported.
# Note: 
#   Not guaranteed to work universally, depends on the specific job.
#   Always check if what this script does makes sense for your purposes.


import os
import sys
import argparse
import glob
import jobCheck

 
if __name__=='__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Job resubmission.')
    parser.add_argument('--jd', required=True,
                        help='Condor job description file managing the submission.')
    parser.add_argument('--dir', default=os.getcwd(),
                        help='Directory to scan for files (default: cwd)')
    parser.add_argument('--starting_tag', default='###starting###',
                        help='Starting tag, default is "###starting###".')
    parser.add_argument('--done_tag', default='###done###',
                        help='Done tag, default is "###done###".')
    parser.add_argument('--ntags', default=None,
                         help='Number of expected starting and done tags per job.')
    args = parser.parse_args()

    # print arguments
    print('running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # some more parsing
    jobfile = os.path.abspath(args.jd)
    if args.ntags is not None: args.ntags = int(args.ntags)    

    # find error log files
    print('finding error log files...')
    condorpattern = os.path.join(args.dir,'*_err_*')
    elfiles = glob.glob(condorpattern)
    nelfiles = len(elfiles)
    print('found {} error log files.'.format(nelfiles))
    print('start scanning...')

    # loop over all error log files found above
    # and find those corresponding to unfinished/failed jobs
    errorfiles = []
    nerrortag = 0
    nunfinished = 0
    for elfile in elfiles:
        if jobCheck.check_start_done(elfile,
                starting_tag = args.starting_tag,
                done_tag = args.done_tag,
                ntarget = args.ntags): 
            errorfiles.append(elfile)
            nunfinished += 1
        elif jobCheck.check_error_content(elfile):
            errorfiles.append(elfile)
            nerrortag += 1
    nerrorfiles = len(errorfiles)
    print('found {} error log files'.format(nunfinished)
            +' corresponding to unfinished jobs and {}'.format(nerrortag)
            +' corresponding to known error content.')

    # handle case of no failures
    if( nerrorfiles==0 ):
        print('nothing to resubmit; exiting.')
        sys.exit()

    # sort the list of error log files
    errorfiles = sorted(errorfiles)
    
    # get corresponding job ids
    jobids = []
    for errorfile in errorfiles:
        jobid = os.path.basename(errorfile).split('_')[3]
        try:
            jobidtest = int(jobid)
            jobids.append(jobid)
        except:
            msg = 'WARNING in jobResubmission.py:'
            msg += ' could not extract job id number'
            msg += ' for the error log file {}.'.format(errorfile)
            msg += ' Does its name have the expected format?'
            print(msg)
            jobids.append(None)

    # check corresponding output log files and find the executables
    shfiles = []
    shtag = '###exename###:'
    for errorfile in errorfiles:
        outputfile = errorfile.replace('_err_', '_out_')
        with open(outputfile,'r') as f:
            line = f.readline()
            if( shtag in line):
                shfile = line.replace(shtag,'').strip(' \t\n')
                shfiles.append(shfile)
            else:
                msg = 'WARNING in jobResubmission.py:'
                msg += ' could not retrieve corresponding executable'
                msg += ' for the output log file {}.'.format(outputfile)
                msg += ' Check if the file contains the tag "{}".'.format(shtag)
                print(msg)
                shfiles.append(None)

    # do some printing and checks
    print('found following resubmission strategy:')
    for i in range(nerrorfiles):
        print(' - file {} (job {}) -> exe {}'.format(errorfiles[i],jobids[i],shfiles[i]))
    print('continue with resubmission?')
    go = raw_input()
    if go!='y': sys.exit()

    # check for earlier warnings and stop if found
    # (later: maybe allow to override this in some cases)
    if( None in jobids or None in shfiles ):
        msg = 'ERROR: cannot proceed with resubmission,'
        msg += ' since not all jobids and/or shfiles'
        msg += ' could be properly retrieved;'
        msg += ' exiting.'
        sys.exit()

    # remove the corresponding jobs
    print('canceling jobs that will be resubmitted...')
    print('(if not done so before)')
    for jobid in jobids:
        os.system('condor_rm {}'.format(jobid))

    # modify the job description file and submit jobs
    print('resubmitting...')
    with open(jobfile,'r') as f:
        lines = f.readlines()
    for shfile in shfiles:
        lines[0] = 'executable = {}\n'.format(shfile)
        with open(jobfile,'w') as f:
            for line in lines: f.write(line)
        os.system('condor_submit {}'.format(jobfile))

    # remove old version of error, output and log files
    print('removing old log files...')
    for i in range(nerrorfiles):
        errorfile = errorfiles[i]
        outputfile = errorfile.replace('_err_', '_out_')
        logfile = errorfile.replace('_err_', '_log_')
        os.system('rm {}'.format(errorfile))
        os.system('rm {}'.format(outputfile))
        os.system('rm {}'.format(logfile))
    
    print('done')
