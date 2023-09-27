###############################################
# select relevant histograms from a ROOT file #
###############################################

import sys
import os
sys.path.append(os.path.abspath('../../tools'))
import histtools as ht
import listtools as lt


def select_histnames(
    inputfile,
    processes=['all'],
    regions=[],
    variablenames=[], 
    includesystematics=None, 
    excludesystematics=None, 
    verbose=False):
  ### get all relevant histogram names from a ROOT file with histograms
  # note: very specific to naming conventions!
 
  # initializations
  doallprocesses = (len(processes)==1 and processes[0]=='all')
  
  # requirement: the histogram name must contain at least one includesystematic (or nominal)
  # and may not contain any of the excludesystematics
  mustcontainone = []
  if includesystematics is not None: mustcontainone = includesystematics + ['nominal']
  maynotcontainone = []
  if excludesystematics is not None: maynotcontainone = excludesystematics
  
  # shortcut requirements for when only one process or variable is requested
  mustcontainall = []
  if( len(processes)==1 and not doallprocesses ): mustcontainall.append(processes[0])
  if( len(variablenames)==1 ): mustcontainall.append(variablenames[0])
  
  # do loading and initial selection
  histnames = ht.loadhistnames(inputfile,
    mustcontainone=mustcontainone,
    maynotcontainone=maynotcontainone,
    mustcontainall=mustcontainall)
  
  # printouts for testing
  if verbose:
    print('Initial selection:')
    print(' - mustcontainone: {}'.format(mustcontainone))
    print(' - mustontainall: {}'.format(mustcontainall))
    print(' - maynotcontainone: {}'.format(excludesystematics))
    print('Resulting number of histograms: {}'.format(len(histnames)))
  
  # select processes
  if not doallprocesses:
    mustcontainone = ['{}_'.format(p) for p in processes]
    histnames = lt.subselect_strings(histnames, mustcontainone=mustcontainone)[1]

  # select processes (special case for processses split at particle level)
  # (to finish implementation later)
  '''if splitvariable is not None and splitprocess is not None:
    mustcontainone = []
    mustcontainone.append( '{}0_'.format(splitprocess) )
    for i in [1,2,3,4,5]:
      mustcontainone.append( '{}{}{}'.format(splitprocess, i, splitvariable.replace("_","")) )
    maynotcontainone=[splitprocess]
    histnames_sig = lt.subselect_strings(histnames, mustcontainone=mustcontainone)[1]
    histnames_back = lt.subselect_strings(histnames, maynotcontainone=maynotcontainone)[1]
    histnames = histnames_sig + histnames_back'''
  
  # select regions  
  mustcontainone = ['_{}_'.format(region) for region in regions]
  histnames = lt.subselect_strings(histnames, mustcontainone=mustcontainone)[1]

  # select variables
  mustcontainone = ['_{}_'.format(variablename) for variablename in variablenames]
  histnames = lt.subselect_strings(histnames, mustcontainone=mustcontainone)[1]

  # printouts for testing
  if verbose:
    print('Further selection (processes, regions and variables):')
    print('Resulting number of histograms: {}'.format(len(histnames)))

  return histnames
