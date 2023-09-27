#!/usr/bin/env python3
# coding: utf-8

from __future__ import annotations

#########################################
# Merge output of eventloop.py for data #
#########################################
# Only needed if primary datasets were processed separately rather than merged.
# In this case, there are still duplicate events in the output files,
# that need to be filtered out.

# Adapted from here (for NanoAOD files):
# https://github.com/GhentAnalysis/nanoSkimming/blob/main/merging/haddnanodata.py

import sys
import os
import math
from functools import partial
from typing import Any
from typing import List
import numpy as np
import awkward as ak
import uproot

try:
    import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def find_tree_names(
    input_paths: List[str]
    ) -> List[str]:
    ### find trees in input files
    # this part depends on the naming convention of trees in the input files.
    treenames = None
    for input_path in input_paths:
        with uproot.open(input_path) as f:
            keys = [key.split(';')[0] for key in f.keys()]
            thistreenames = sorted([key for key in keys if key.endswith('Events')])
        if treenames is None: treenames = thistreenames[:]
        elif treenames!=thistreenames:
            msg = 'ERROR: incompatible sets of tree names:\n'
            msg += '{}\n{}'.format(treenames, thistreenames)
    # temporary to avoid errors
    #treenames = [t for t in treenames if 'cfcontrolregion' not in t]
    return treenames


def hadddata(
    output_path: str,
    input_paths: List[str],
    force: bool = False,
    keep_branches: List[str] | None = None,
    step_size: int = 100000,
    verbose: bool = False,
) -> None:

    # expand variables
    expand = lambda path: os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
    input_paths = [expand(input_path) for input_path in input_paths]
    output_path = expand(output_path)

    # check if all input files exist
    for input_path in input_paths:
        if not os.path.isfile(input_path):
            msg = 'ERROR: input file {} does not seem to exist.'.format(input_path)
            raise Exception(msg)

    # define index columns and check if they are contained in branches to keep
    index_columns = ["event", "run", "luminosityBlock"]
    if keep_branches is not None:
        for index_column in index_columns:
            if index_column not in keep_branches:
                msg = 'ERROR: keep_branches does not contain required index branches'
                msg += ' {}.'.format(index_columns)
                raise Exception(msg)

    # prepare the output
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    elif os.path.exists(output_path):
        if not force:
            msg = 'ERROR: output file {} already exists.'.format(output_path)
            msg += ' Either remove it or run with force=True to overwrite it.'
            raise Exception(msg)
        else:
            msg = 'WARNING: overwriting existing file {}...'.format(output_path)
            print(msg)
            os.remove(output_path)
    output_file = uproot.create(output_path)

    # loop over tree names
    treenames = sorted(find_tree_names(input_paths))
    if verbose:
        print('Found following tree names in input files:')
        for treename in treenames: print('  - {}'.format(treenames))
    for treename in treenames:
        if verbose: print('Now running on tree {}'.format(treename))

        # get trees from input files
        # note: it is convenient to make sure that the first tree in the list is not empty
        #       in order to correctly initialize the output tree and safely use extend.
        # note: remaining problem: if all trees are empty, the output tree is not written.
        trees = [uproot.open(input_path)[treename] for input_path in input_paths]
        temp = [tree.num_entries for tree in trees]
        sorted_ids = np.argsort(temp)[::-1]
        trees = [trees[idx] for idx in sorted_ids]
        tree1 = trees[0]

        # read index columns over the full reference file
        index = tree1.arrays(index_columns)

        # prepare counts
        n_written = 0
        n_overlap = 0

        # iteration helper
        def iterate(tree, ntree, ntrees):
            if verbose:
                print(f"  Iterating through tree {ntree}/{ntrees} with {tree.num_entries} events")
            # safety for empty trees
            if tree.num_entries==0:
                #chunk = {}
                #for key in tree.keys(): chunk[key] = np.zeros(0)
                #chunk = ak.Array(chunk)
                #return [chunk]
                return []
            progress = (
                partial(tqdm.tqdm, total=int(math.ceil(tree.num_entries / step_size)))
                if( verbose and HAS_TQDM ) else (lambda gen: gen) )
            return progress(tree.iterate(step_size=step_size, filter_name=keep_branches))

        # fill chunks of the first tree
        for chunk in iterate(tree1, 1, len(trees)):
            # update counts
            n_written += len(chunk)
            # extend the output tree
            chunk = dict(zip(chunk.fields, ak.unzip(chunk)))
            output_keys = [key.split(';')[0] for key in output_file.keys()]
            if treename in output_keys: output_file[treename].extend(chunk)
            else: output_file[treename] = chunk

        # fill chunks of the other trees
        for idx, tree in enumerate(trees[1:]):
            for chunk in iterate(tree, idx+2, len(trees)):
                # determine a mask of events in tree that are also in tree1
                mask = np.isin(chunk[index_columns], index, assume_unique=True)
                chunk = chunk[~mask]
                # update counts
                n_written += len(chunk)
                n_overlap += ak.sum(mask)
                # skip the chunk if all events are overlapping
                if ak.all(mask): continue
                # extend the output tree
                chunk = dict(zip(chunk.fields, ak.unzip(chunk)))
                output_file[treename].extend(chunk)
                # update the index
                chunkindex = {key: chunk[key] for key in index_columns}
                chunkindex = ak.Array(chunkindex)
                index = ak.concatenate((index, chunkindex))

        if verbose:
            print(f"  Written {n_written} and found {n_overlap} overlapping event(s)")


if __name__ == "__main__":
   
    # read command line arguments 
    import argparse
    parser = argparse.ArgumentParser(
        description="joins data files and removes duplicate events")
    parser.add_argument("-o", "--outputfile",
        help="Path to the output file to be created")
    parser.add_argument("-i", "--inputfiles", nargs='+',
        help="Path to input files to merge")
    parser.add_argument("-f", "--force", default=False, action='store_true',
        help="Whether to overwrite output file if it already exists")
    parser.add_argument("--step-size", "-s", type=int, default=100000,
        help="step size for iterations; default: 100000")
    parser.add_argument("--verbose", "-v", default=False, action="store_true",
        help="verbose output, potentially with tqdm if installed")
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg))) 

    keep_branches = None

    # merge the files
    hadddata(
        output_path=args.outputfile,
        input_paths=args.inputfiles,
        force=args.force,
        keep_branches=keep_branches,
        step_size=args.step_size,
        verbose=args.verbose )
