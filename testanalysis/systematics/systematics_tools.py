#################################################
# Tools for formatting / processing systematics #
#################################################

import sys
import os
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
import objectselection.jetuncs as jetuncs
from systematics_type import systematics_type


def get_selection_systematics(events, systematics, includenominal=False):
    ### get a dict of selection systematics
    # input arguments:
    # - systematics: a list of valid systematic names
    #   (must be keys in the systematics_type dict)
    # - includenominal: include "nominal" in output
    # returns:
    # a dict with all selection systematics mapped to their variations,
    # e.g. ["jec"] -> {"jec": ["jesTotalUp", "jesTotalDown", ... ]}.
    # if there are no variations, the key is repeated as value,
    # e.g. {"nominal": ["nominal"]}
    selection_systematics = ([systematic for systematic in systematics
      if systematics_type[systematic].endswith('selection')])
    selection_systematics = {key: [] for key in selection_systematics}
    for selection_systematic in selection_systematics.keys():
        if selection_systematic == 'jec':
            variations = jetuncs.get_available_jec_variations(events)
            for variation in variations:
                selection_systematics[selection_systematic].append(variation+'Up')
                selection_systematics[selection_systematic].append(variation+'Down')
        elif selection_systematic == 'jer':
            variations = jetuncs.get_available_jer_variations(events)
            for variation in variations:
                selection_systematics[selection_systematic].append(variation+'Up')
                selection_systematics[selection_systematic].append(variation+'Down')
        elif selection_systematic == 'uncl':
            selection_systematics[selection_systematic].append('unclustEnUp')
            selection_systematics[selection_systematic].append('unclustEnDown')
        else:
            msg = 'ERROR: unrecognized selection systematic: {}'.format(selection_systematic)
            raise Exception(msg)
    if includenominal: selection_systematics['nominal'] = ['nominal']
    return selection_systematics

def get_weight_systematics(reweighter, systematics, includenominal=False):
    ### get a dict of weight systematics
    # input arguments:
    # - reweighter: a CombinedReweighter instance
    # - systematics: a list of valid systematic names
    #   (must be keys in the systematics_type dict)
    # - includenominal: include "nominal" in output
    # returns:
    # a dict with all weight systematics mapped to their variations,
    # e.g. ["muonid"] -> {"muonid": ["stat_up", "stat_down", "syst_up", "syst_down"]}.
    # if there are no variations, the key is repeated as value,
    # e.g. {"nominal": ["nominal"]}
    weight_systematics = ([systematic for systematic in args.systematics
      if systematics_type[systematic]=='weight'])
    weight_systematics = {key: [] for key in weight_systematics}
    for weight_systematic in weight_systematics.keys():
        weight_systematics[weight_systematic] = reweighter.get_variations(weight_systematic)
    return weight_systematics
