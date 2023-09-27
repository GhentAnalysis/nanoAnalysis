########################
# Test histogram tools #
########################

import sys
import os
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
import tools.histtools2 as ht
from tools.histogram import Histogram


if __name__=='__main__':

    histfile = sys.argv[1]

    # load all histogram names
    histnames = ht.loadallhistnames(histfile)
    print(len(histnames))
    print(histnames[:5])

    # load selected histogram names
    histnames = ht.loadhistnames(histfile, mustcontainall=['nominal'])
    print(len(histnames))
    print(histnames[:5])

    # load histograms
    hists = ht.loadhistogramlist(histfile, histnames)
    print(len(hists))
    hist = hists[0]
    print(hist)
    print(type(hist))
    hist = Histogram.from_uproot(hist)
    print(hist)
    print(type(hist))

    # change histogram name
    print(hist.name)
    hist.name = 'test'
    print(hist.name)

    # change histogram values
    print(hist.values)
    hist.values[1] = 10.
    print(hist.values)
