######################################################################
# Tools for event selection using trigger and MET filter information #
######################################################################
# note: largely copied/translated from 
# https://github.com/LukaLambrecht/ewkino/blob/ttW/ttWAnalysis/eventselection/src/eventSelections.cc

# imports
import awkward as ak
import numpy as np

# external functions
# (meant to be called from outside while doing event selections)

def pass_met_filters(events, filters=None):
    if filters is None: filters = ['METFilters']
    mask = (events.event >= 0)
    for f in filters:
        # check if filter is present in events
        if f not in events.Flag.fields:
            msg = 'ERROR: filter {} not found.'.format(f)
            raise Exception(msg)
        mask = (mask & events.Flag[f])
    return mask

def pass_any_lepton_trigger(events, triggers=None):
    if triggers is None:
        triggers = [
          'trigger_e',
          'trigger_ee',
          'trigger_eee',
          'trigger_m',
          'trigger_mm',
          'trigger_mmm',
          'trigger_em',
          'trigger_eem',
          'trigger_emm',
          'trigger_emm'
        ]
    mask = (events.event < 0)
    for trigger in triggers:
        # check if trigger is present in events
        if trigger not in events.HLT.fields:
            msg = 'ERROR: trigger {} not found.'.format(trigger)
            raise Exception(msg)
        mask = (mask | events.HLT[trigger])
    return mask
