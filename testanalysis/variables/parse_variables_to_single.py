###################################################
# parse a json list of double variables to single #
###################################################
# use case:
# - make modifications to main list of double variables,
#   with event BDT as primary variable, and other as secondary
#   (e.g. add variables, change binning, etc).
# - use this script to read the double variables
#   and make a list of all secondary as single variables
#   (e.g. to run the response matrices).
# - this will make sure they are exactly consistent
#   without manual tuning and copying.

# remark: run with python3 instead of python2
# to avoid problems with unicode vs string!

import sys
import os
sys.path.append('../../Tools/python')
import variabletools as vt

if __name__=='__main__':

  inputfile = sys.argv[1]
  outputfile = sys.argv[2]
  
  varlist = vt.read_variables( inputfile )
  newvarlist = []
  for var in varlist:
    newvarlist.append( var.secondary )
  vt.write_variables_json( newvarlist, outputfile, builtin=False )
