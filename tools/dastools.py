####################################
# Tools for reading files from DAS #
####################################

import sys
import os


def get_sample_files( datasetname, 
                      filemode='das',
                      privateprod=False,
                      redirector='root://cms-xrd-global.cern.ch/',
                      istest=False,
                      maxfiles=None,
		      verbose=False ):
  ### get a list of input files from a sample name
  # input arguments:
  # - datasetname: name of the data set on DAS (for filemode 'das')
  #   OR name of the folder holding input files (for filemode 'local')
  #   OR str containing a comma-separated list of file names (on DAS or locally according to filemode))
  #   (note: interpreted as list of file names if a comma is present, directory or dataset otherwise!)
  # - filemode: choose from 'das' or 'local';
  #   in case of 'das', will read all files belonging to the specified dataset from DAS;
  #   in case of 'local', will read all files in the specified folder on the local filesystem.)
  # - privateprod: whether the dataset on DAS was privately produced (ignored in filemode 'local')
  # - redirector: redirector used to access remote files (ignored in filemode 'local'))
  # - istest: return only first file (for testing)
  # - maxfiles: return only specified number of first files
  # note: the DAS client requires a valid proxy to run,
  #       set it before calling this function with set_proxy() (see below)

  # check if a directory is provided or a list of filenames
  runfilesearch = True
  if ',' in datasetname: runfilesearch = False

  # parse the provided redirector
  redirector = redirector.rstrip('/')+'/'

  # make a list of input files
  if runfilesearch:
    # make a list of input files based on provided directory or dataset name,
    # details depend on the chosen filemode
    if filemode=='das':
      # make and execute the DAS client command
      if verbose: print('running DAS client to find files in dataset {}...'.format(datasetname))
      instance = ''
      if privateprod: instance = ' instance=prod/phys03'
      dascmd = "dasgoclient -query 'file dataset={}{}' --limit 0".format(datasetname,instance)
      dasstdout = os.popen(dascmd).read()
      dasfiles = sorted([el.strip(' \t') for el in dasstdout.strip('\n').split('\n')])
      if verbose:
        print('DAS client ready; found following files ({}):'.format(len(dasfiles)))
        for f in dasfiles: print('  - {}'.format(f))
      inputfiles = [redirector+f for f in dasfiles]
    elif filemode=='local':
      # read all root files in the given directory
      inputfiles = ([os.path.join(datasetname,f) for f in os.listdir(datasetname)
                     if f[-5:]=='.root'])
  else:
    # parse the provided comma-separated list into a list
    inputfiles = [el for el in datasetname.split(',') if len(el)!=0]
    if( filemode=='das' and redirector is not None ):
      for i,inputfile in enumerate(inputfiles):
        # check if the file name has a redirector already
        if 'root://cms-xrd' in inputfile: continue
        # add the redirector
        inputfiles[i] = redirector+inputfile

  # check number of input files
  if len(inputfiles)==0:
    raise Exception('ERROR: list of input files is empty.')
  if istest:
    print('WARNING: running in test mode, only one file will be processed.')
    inputfiles = [inputfiles[0]]
  if( maxfiles is not None and maxfiles>0 and maxfiles<len(inputfiles) ):
    print('WARNING: returning only {} out of {} files.'.format(maxfiles,len(inputfiles)))
    inputfiles = inputfiles[:maxfiles]

  return inputfiles


def export_proxy( proxy ):
  ### export a provided proxy to the system variables
  print('exporting proxy to {}'.format(proxy))
  os.environ["X509_USER_PROXY"] = proxy
