#######################
# convert json to txt #
#######################

import sys
import os
sys.path.append('../../Tools/python')
import variabletools as vt

if __name__=='__main__':

  inputfile = sys.argv[1]
  outputfile = inputfile.replace('.json','.txt')
  
  if os.path.exists(outputfile):
    raise Exception('ERROR: output file already exists.')

  varlist = vt.read_variables( inputfile )
  vt.write_variables_txt( varlist, outputfile ) 
