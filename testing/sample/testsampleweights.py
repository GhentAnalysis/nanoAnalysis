############################################
# Test reading of sample generator weights #
############################################

import sys
import os
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
from samples.sampleweights import SampleWeights

if __name__=='__main__':

    # read input file from command line
    inputfile = sys.argv[1]
    
    # make a SampleWeights object
    weights = SampleWeights(inputfile)
    print(weights)

    # print pdf and scale variations
    print(weights.pdfweights(returntype='minmax'))
    print(weights.scaleweights(returntype='all'))
    print(weights.scaleweights(returntype='minmax'))
    print(weights.scaleweights(returntype='envelope'))
