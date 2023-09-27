##########################################
# Add custom variables to lepton objects #
##########################################

# imports
import sys
import os
import awkward as ak

def add_electron_variables(events, variables=['all']):
    ### add electron variables to events
    # input arguments:
    # - events: an object of type NanoEventsArray
    # - variables: list of variable names to add
    #   (default: all variables defined in this function)
    allvariables = ('all' in variables)
    checkvars = variables[:]
    # ptRatio
    # defined as 1/(relative isolation in matched jet + 1).
    # note: if no matched jet, Electron.jetRelIso is set to pfRelIso04_all,
    #       so it always has a defined value (i.e. not None).
    if( allvariables or 'jetPtRatio' in variables):
        if not allvariables: checkvars.remove('jetPtRatio')
        events.Electron = ak.with_field(events.Electron,
          1./(events.Electron.jetRelIso+1),
          where='jetPtRatio')
    # deep flavor score of closest jet
    # note: if no matched jet, this value is set to 0 by hand;
    #       this appears to be needed to be able to write the events to a file
    #       (gives errors if this variable contains None).
    if( allvariables or 'jetBTagDeepFlavor' in variables):
        if not allvariables: checkvars.remove('jetBTagDeepFlavor')
        events.Electron = ak.with_field(events.Electron,
          ak.fill_none(events.Electron.matched_jet.btagDeepFlavB,0),
          where='jetBTagDeepFlavor')
    # also update key access
    # note: appears to be needed to fully propagate the addition of this variable.
    events['Electron'] = events.Electron
    # check unrecognized variables
    if( not allvariables and len(checkvars)>0 ):
        msg = 'ERROR in leptonvariables.py / add_electron_variables:'
        msg += ' following variables were not recognized: {}.'.format(checkvars)
        raise Exception(msg)

def add_muon_variables(events, variables=['all']):
    ### add muon variables to events
    # input arguments:
    # - events: an object of type NanoEventsArray
    # - variables: list of variable names to add
    #   (default: all variables defined in this function)
    allvariables = ('all' in variables)
    checkvars = variables[:]
    # ptRatio
    # defined as 1/(relative isolation in matched jet + 1).
    # note: if no matched jet, Electron.jetRelIso is set to pfRelIso04_all,
    #       so it always has a defined value (i.e. not None).
    if( allvariables or 'jetPtRatio' in variables):
        if not allvariables: checkvars.remove('jetPtRatio')
        events.Muon = ak.with_field(events.Muon,
          1./(events.Muon.jetRelIso+1),
          where='jetPtRatio')
    # deep flavor score of closest jet
    # note: if no matched jet, this value is set to 0 by hand;
    #       this appears to be needed to be able to write the events to a file
    #       (gives errors if this variable contains None).
    if( allvariables or 'jetBTagDeepFlavor' in variables):
        if not allvariables: checkvars.remove('jetBTagDeepFlavor')
        events.Muon = ak.with_field(events.Muon,
          ak.fill_none(events.Muon.matched_jet.btagDeepFlavB,0),
          where='jetBTagDeepFlavor')
    # also update key access
    # note: appears to be needed to fully propagate the addition of this variable.
    events['Muon'] = events.Muon
    # check unrecognized variables
    if( not allvariables and len(checkvars)>0 ):
        msg = 'ERROR in leptonvariables.py / add_muon_variables:'
        msg += ' following variables were not recognized: {}.'.format(checkvars)
        raise Exception(msg)
