#############################################
# some small tools for comparing root files #
#############################################
# use case: compare output file before and after code changes,
# when nothing is supposed to have changed.
# so far only works for files containing TH1 (no trees or other objects)

import os
import sys
import argparse
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
import tools.histtools as ht

def th1_equal(hist1, hist2, threshold=1e-12):
    ### compare bin contents of two TH1 objects
    nbins = hist1.GetNbinsX()
    if( hist2.GetNbinsX()!=nbins ): return False
    for i in range(0,nbins+2):
        if( abs(hist1.GetBinContent(i)-hist2.GetBinContent(i))>threshold ): return False
        if( abs(hist1.GetBinError(i)-hist2.GetBinError(i))>threshold ): return False
    return True


if __name__=='__main__':

    # parse arguments
    parser = argparse.ArgumentParser('Compare two root files containing TH1')
    parser.add_argument('--file1', required=True, type=os.path.abspath)
    parser.add_argument('--file2', required=True, type=os.path.abspath)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
      print('  - {}: {}'.format(arg,getattr(args,arg)))

    # compare histogram names
    print('Checking histogram names...')
    hnames1 = sorted( ht.loadallhistnames(args.file1) )
    hnames2 = sorted( ht.loadallhistnames(args.file2) )
    print('Found {} histograms in file 1 and {} in file 2.'.format(len(hnames1),len(hnames2)))
    if hnames1!=hnames2:
        print('Lists of histogram names are different!')
        sys.exit()
    print('Lists of histogram names are equal.')

    # compare histogram contents
    print('Checking histogram bin contents...')
    hists1 = ht.loadhistogramlist(args.file1, hnames1)
    hists2 = ht.loadhistogramlist(args.file2, hnames1)
    diff = False
    for hist1, hist2 in zip(hists1, hists2):
        if not th1_equal(hist1, hist2):
            diff = True
            print('Difference found in histogram {}!'.format(hist1.GetName()))
    if not diff:
        print('Bin contents are equal.')
