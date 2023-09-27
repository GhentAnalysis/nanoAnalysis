########################################
# test the reading of a variables json #
########################################

import sys
import os
sys.path.append('../../Tools/python')
import variabletools as vt

if __name__=='__main__':

  inputfiles = sys.argv[1:]

  for inputfile in inputfiles:
    print('Reading variable list {}...'.format(inputfile))
    varlist = vt.read_variables( inputfile )
    for var in varlist:
      print(var)
