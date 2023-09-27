##########################################
# Get yields of histograms in ROOT files #
##########################################

import sys
import os
import argparse
sys.path.append(os.path.abspath('../../tools'))
import histtools as ht


if __name__=='__main__':

    # parse arguments
    parser = argparse.ArgumentParser(description='Print yields')
    parser.add_argument('-i', '--inputfiles', required=True, type=os.path.abspath, nargs='+')
    parser.add_argument('-n', '--histname', required=True)
    parser.add_argument('-s', '--sorting', default=None,
      choices=[None, 'ascending', 'descending', 'alpha'])
    args = parser.parse_args()

    # print arguments
    #print('Running with following configuration:')
    #for arg in vars(args):
    #    print('  - {}: {}'.format(arg,getattr(args,arg)))

    # initializations
    yields = {}

    # loop over input files
    for inputfile in args.inputfiles:
        # get the histogram
        hist = ht.loadhistogramlist(inputfile, [args.histname])
        if len(hist)!=1: continue
        hist = hist[0]
        integral = hist.Integral()
        yields[inputfile] = integral

    # sort the results
    if args.sorting == 'alpha': yields = sorted(yields.items())
    elif args.sorting == 'ascending':
        yields = sorted(yields.items(), key=lambda x:x[1])
    elif args.sorting == 'descending':
        yields = sorted(yields.items(), key=lambda x:x[1], reverse=True)
    else: yields = yields.items()

    # print results
    for f, y in yields:
        print(' - {}: {}'.format(f, y))
