################################
# Add combinations of triggers #
################################

# imports
import os
import sys
import json
import awkward as ak

def load_triggerdefs():
    ### internal helper function to load the json file with trigger definitions
    # path to json file is hard-coded for now, maybe later use an argument.
    triggerfile = os.path.join(os.path.dirname(__file__), 'triggers.json')
    with open(triggerfile) as f:
        triggerdefs = json.load(f)
    return triggerdefs

def add_trigger_variables(events, year=None, returntriggerdefs=False):
    ### add combined trigger definitions to NanoEvents
    # load trigger definitions
    if year is None:
        raise Exception('ERROR in add_trigger_variables: year is not set.')
    triggerdefs = load_triggerdefs()
    if year not in triggerdefs.keys():
        raise Exception('ERROR in add_trigger_variables: year {} not recognized;'.format(year)
          +' options are: {}'.format(triggerdefs.keys()))
    triggerdefs = triggerdefs[year]
    # loop over combined triggers
    for trigger, hlts in triggerdefs.items():
        # make a mask for this combined trigger
        mask = (events.event < 0)
        for hlt in hlts:
            # check if trigger is present in events
            if hlt not in events.HLT.fields:
                msg = 'WARNING: trigger {} not found'.format(hlt)
                msg += ' (needed for definition of {});'.format(trigger)
                msg += ' this trigger will be ignored in the combination.'
                print(msg)
                continue
            mask = (mask | events.HLT[hlt])
        # add this as a new variable
        events.HLT = ak.with_field(events.HLT, mask, where=trigger)
        events['HLT'] = events.HLT
    # return added triggers if requested
    if returntriggerdefs: return triggerdefs
