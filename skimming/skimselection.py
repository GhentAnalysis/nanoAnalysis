#################################################
# Collection of selection functions for skimmer #
#################################################
# Use by importing the function skimselection (see bottom).

import sys
import os
import awkward as ak

def nlightlepton(events, n,
                 selected_muons=None,
                 selected_electrons=None):
    ### perform >= n light lepton selection
    if selected_muons is None: selected_muons = events.Muon
    if selected_electrons is None: selected_electrons = events.Electron
    selection = (ak.num(selected_muons.pt) + ak.num(selected_electrons.pt) >= n)
    return selection

def dilightlepton(events, **kwargs):
    ### perform >= 2 light lepton selection
    return nlightlepton(events, 2, **kwargs)

def trilightlepton(events, **kwargs):
    ### perform >= 3 light lepton selection
    return nlightlepton(events, 3, **kwargs)

def multilightlepton(events,
                     selected_muons=None,
                     selected_electrons=None):
    ### perform selection on either > 2 light leptons or 2 same-sign light leptons
    if selected_muons is None: selected_muons = events.Muon
    if selected_electrons is None: selected_electrons = events.Electron
    nlep = ak.num(selected_muons.pt) + ak.num(selected_electrons.pt)
    lepcharge = ak.concatenate((selected_muons.charge, selected_electrons.charge), axis=1)
    lepchargesum = abs(ak.sum(lepcharge, axis=1))
    selection = ( (nlep > 2) | ((nlep == 2) & (lepchargesum == 2)) )
    return selection

def skimselection(events, 
                  selectionid=None,
                  muon_mask=None,
                  electron_mask=None):
    ### perform event selection for skimming
    # input arguments:
    # - events: object of type NanoEventsArray
    # - selectionid: selection identifier string
    # - muon_mask: mask to define good muons
    # - electron_mask: mask to define good electrons

    # do object selection
    selected_muons = events.Muon
    if muon_mask is not None: selected_muons = events.Muon[muon_mask]
    selected_electrons = events.Electron
    if electron_mask is not None: selected_electrons = events.Electron[electron_mask]

    # do event selection
    if( selectionid=='dilightlepton' ): 
        return dilightlepton(events,
                 selected_muons=selected_muons, 
                 selected_electrons=selected_electrons)
    if( selectionid=='trilightlepton' ):
        return trilightlepton(events,
                 selected_muons=selected_muons,
                 selected_electrons=selected_electrons)
    if( selectionid=='multilightlepton' ):
        return multilightlepton(events,
                 selected_muons=selected_muons,
                 selected_electrons=selected_electrons)

    # raise error if selection parameters are invalid
    msg = 'ERROR in skimselection:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)
