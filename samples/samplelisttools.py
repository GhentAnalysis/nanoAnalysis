##################################
# Tools for reading sample lists #
##################################

import sys
import os
import argparse
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
import tools.argparsetools as apt
from samples.sample import Sample, SampleCollection


def readsamplelist( samplelistpaths, sampledir=None, doyear=True ):
    ### returns a SampleCollection from a list of sample list files
    # input arguments:
    # - samplelistpaths: either string or list of strings 
    #                    representing path(s) to sample list file(s).
    # - sampledir: path to directory containing the samples.
    #              if not specified, the path attribute of each sample is None
    #              and no check on sample existence is done.
    # - doyear: extract sample year from file name.
    # note: each line of the sample list is assumed to be of the form
    #       <process_name> <sample file> <cross_section>
    #       where the individual elements are separated by spaces 
    #       and should not contain spaces,
    #       and the cross_section is optional (defaults to 0);
    
    collection = SampleCollection()
    if isinstance(samplelistpaths, str): samplelistpaths = [samplelistpaths]

    for s in samplelistpaths:
        if not os.path.exists(s):
            raise Exception('ERROR in readsamplelist:'
                    +' sample list {} does not exist.'.format(s))
    if( sampledir is not None and not os.path.exists(sampledir) ):
        raise Exception('ERROR in readsamplelist:'
                    +' sample directory {} does not exist.'.format(sampledir))

    collection.read_from_files( samplelistpaths, 
                sampledir=sampledir,
                doyear=doyear,
                suppress_exception=True )

    # check non-existing samples
    if sampledir is not None:
        missing_samples = collection.check_paths( return_list=True )
        if len(missing_samples)>0:
            msg = 'ERROR in readsamplelist:'
            msg += ' the following samples were not found in {}:\n'.format(sampledir)
            for s in missing_samples: msg += '{}\n'.format(s)
            raise Exception(msg)

    # check inconsistent years
    if doyear:
        years = set([s.year for s in collection.get_samples()])
        if len(years)>1:
            msg = 'ERROR in readsamplelist:'
            msg += ' found multiple years in sample list: {}'.format(years)
            raise Exception(msg)

    return collection


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Read sample list')
    parser.add_argument('-l', '--samplelist', required=True, type=os.path.abspath)
    parser.add_argument('-d', '--sampledir', type=apt.path_or_none, default=None)
    args = parser.parse_args()

    samples = readsamplelist( args.samplelist, sampledir=args.sampledir )
    print(samples)
