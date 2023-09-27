##############################
# Run full chain of commands #
##############################
# How to use:
# - check that all options in the separate steps are set correctly.
# - run with "nohup python3 chain.py &> [name of a log file] &".
#   some explanation:
#   - "python3 chain.py" trivially runs this script.
#   - "&> [name of a log file]" redirects both stdout and stderr to the log file you chose.
#   - "nohup [...] &" runs the command in the background, 
#     and keeps it running even after you close the terminal.
# - alternatively, you could open a screen session and simply run 
#   "python3 chain.py &> [name of a log file]".
#   this gives you more control over terminating the command,
#   since with "nohup &", the process ID seems to be irretrievable 
#   after logging out of the m-machine (?).


import sys
import os
import time
sys.path.append(os.path.abspath('../jobsubmission'))
import jobtools


def wait(**kwargs):
    while True:
        if not jobtools.jobs_are_running(verbose=True, **kwargs): return
        time.sleep(60)
    
# set output directory
outputdir = sys.argv[1]

# run event loop and wait for jobs to finish
cmd = 'python3 eventloop_loop.py'
cmd += ' {}'.format(outputdir)
os.system(cmd)
wait()

# run data merging
cmd = 'python3 mergedatatrees_loop.py'
cmd += ' {}'.format(outputdir)
os.system(cmd)
wait()

# run binning and wait for jobs to finish
cmd = 'python3 binner_loop.py'
cmd += ' {}'.format(outputdir)
os.system(cmd)
wait()

# run merging folders
# (run locally!)
cmd = 'python3 mergedatasim.py'
cmd += ' {}'.format(outputdir)
os.system(cmd)

# run merging histograms
# (run locally!)
cmd = 'python3 mergehists_loop.py'
cmd += ' {}'.format(outputdir)
os.system(cmd)
