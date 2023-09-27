###################################################
# parse a json list of single variables to double #
###################################################
# use case:
# - opposite of parse_variables_to_single!
# - make modifications to main list of single variables,
#   (e.g. add variables, change binning, etc).
# - use this script to read the single variables
#   and make a list of double variables,
#   with the read variables as secondary,
#   and a variable specified below as primary.
# - this will make sure they are exactly consistent
#   without manual tuning and copying.

import sys
import os
sys.path.append('../../Tools/python')
import variabletools as vt
from variabletools import HistogramVariable, DoubleHistogramVariable

if __name__=='__main__':

  inputfile = sys.argv[1]
  outputfile = sys.argv[2]

  primary = ({
    "name": "_eventBDT",
    "variable": "_eventBDT",
    "bins": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    "axtitle": "BDT output score"
  })
  primary = HistogramVariable.fromdict(primary)
 
  varlist = vt.read_variables( inputfile )
  newvarlist = []
  for var in varlist:
    name = primary.name + var.name.strip('_')
    newvarlist.append( DoubleHistogramVariable( name, primary, var ) )
  vt.write_variables_json( newvarlist, outputfile, builtin=False )
