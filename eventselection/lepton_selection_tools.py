######################################################
# Tools for event selection using lepton information #
######################################################
# note: largely copied/translated from 
# https://github.com/LukaLambrecht/ewkino/blob/ttW/ttWAnalysis/eventselection/src/eventSelections.cc

# imports
import awkward as ak
import numpy as np

# internal helper functions
# (just for convenience in the functions below,
#  not meant to be called from outside)

def get_leptons(events,
    electron_mask=None, muon_mask=None):
    ### get selected leptons
    electrons = events.Electron
    if electron_mask is not None: electrons = events.Electron[electron_mask]
    muons = events.Muon
    if muon_mask is not None: muons = events.Muon[muon_mask]
    return (electrons, muons)

# external functions
# (meant to be called from outside while doing event selections)

def clean_leptons_and_jets(events):
    pass # not yet implemented

def pass_lepton_pt_thresholds(events,
    pt_thresholds=None, **kwargs):
    ### check if (selected) leptons pass pt thresholds
    if pt_thresholds is None: return True
    # (note: still to test if the above throws errors,
    #  but in any case it is not expected to be a valid use case.)
    (electrons, muons) = get_leptons(events, **kwargs)
    # get the sorted lepton pts
    pts = ak.sort(ak.concatenate((electrons.pt, muons.pt), axis=1), axis=1, ascending=False)
    # get the leading lepton pts (pad with zeros if not enough leptons in event)
    pts = ak.fill_none(ak.pad_none(pts, len(pt_thresholds), clip=True), 0.)
    # create the mask
    lmask = (np.greater(ak.to_numpy(pts), np.array(pt_thresholds)))
    mask = lmask[:,0] 
    for i in range(1,len(pt_thresholds)): mask = mask & lmask[:,i]
    return mask

def pass_all_leptons_prompt(events, **kwargs):
    ### check if all (selected) leptons in an event are prompt
    (electrons, muons) = get_leptons(events, **kwargs)
    isnonprompt = ak.concatenate((~electrons.isPrompt, ~muons.isPrompt), axis=1)
    mask = (ak.sum(isnonprompt,axis=1)==0)
    return mask

def pass_all_leptons_charge(events, **kwargs):
    ### check if all (selected) leptons are not charge flipped
    (electrons, muons) = get_leptons(events, **kwargs)
    ischargeflip = ak.concatenate((electrons.isChargeFlip, muons.isChargeFlip), axis=1)
    mask = (ak.sum(ischargeflip,axis=1)==0)
    return mask

def get_lepton_invmass(events, **kwargs):
    ### get invariant mass of lepton system
    (electrons, muons) = get_leptons(events, **kwargs)
    leptons = ak.with_name(ak.concatenate((electrons,muons),axis=1), 'PtEtaPhiMCandidate')
    return leptons.sum(axis=1).mass

def get_lepton_invmass_pairs(events, **kwargs):
    ### get invariant masses of all pairs of leptons
    (electrons, muons) = get_leptons(events, **kwargs)
    leptons = ak.with_name(ak.concatenate((electrons,muons),axis=1), 'PtEtaPhiMCandidate')
    pairs = ak.combinations(leptons, 2, fields=['i0', 'i1'])
    return (pairs.i0 + pairs.i1).mass

def pass_mll_lowmass_veto(events, threshold=12,
    electron_mask=None, muon_mask=None):
    ### veto resonances at low invariant mass for same-flavor leptons
    # first for electrons (ignore all muons)
    emasses = get_lepton_invmass_pairs(events,
      electron_mask=electron_mask, muon_mask=(events.Muon.pt < 0.))
    emask = ak.all(emasses >= threshold, axis=1)
    # now for muons (ignore all electrons)
    mmasses = get_lepton_invmass_pairs(events,
      electron_mask=(events.Electron.pt < 0.), muon_mask=muon_mask)
    mmask = ak.all(mmasses >= threshold, axis=1)
    return (emask & mmask)

def pass_tight_lepton_selection(events, nleptons, selectiontype,
    electron_base_mask=None, muon_base_mask=None,
    electron_tight_mask=None, muon_tight_mask=None):
    ### perform tight lepton selection for different selection types
    mask = None
    if(selectiontype=='tight'):
        # normal selection of tight leptons for data vs simulation
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)==nleptons)
    elif(selectiontype=='prompt'):
        # selection of tight prompt leptons (for nonprompt from data)
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)==nleptons)
        if(events.metadata['dtype']=='sim'):
            mask = (mask & pass_all_leptons_prompt(events, 
                electron_mask = electron_tight_mask, muon_mask = muon_tight_mask) )
    elif(selectiontype=='chargegood'):
        # selection of tight leptons with correct charge (for charge flips from data)
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)==nleptons)
        if(events.metadata['dtype']=='sim'):
            mask = (mask & pass_all_leptons_charge(events,
                electron_mask = electron_tight_mask, muon_mask = muon_tight_mask) )
    elif(selectiontype=='irreducible'):
        # combination of prompt and chargegood (for both nonprompt and charge flips from data)
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)==nleptons)
        if(events.metadata['dtype']=='sim'):
            mask = (mask & pass_all_leptons_prompt(events,
                electron_mask = electron_tight_mask, muon_mask = muon_tight_mask) )
            mask = (mask & pass_all_leptons_charge(events,
                electron_mask = electron_tight_mask, muon_mask = muon_tight_mask) )
    elif(selectiontype=='fakerate'):
        # selection of at least one non-tight lepton (for nonprompt from data)
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)<nleptons)
        if(events.metadata['dtype']=='sim'):
            mask = (mask & pass_all_leptons_prompt(events,
                electron_mask = electron_base_mask, muon_mask = muon_base_mask) )
    elif(selectiontype=='efakerate'):
        # same as above but electron fake only (for splitting in muon and electron fakes)
        mask = pass_tight_electron_selection(events, nleptons, 'fakerate',
          electron_base_mask=electron_base_mask, muon_base_mask=muon_base_mask,
          electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
        # to do!
    elif(selectiontype=='mfakerate'):
        # same as above but muon fakes only (for splitting in muon and electron fakes)
        mask = pass_tight_electron_selection(events, nleptons, 'fakerate',
          electron_base_mask=electron_base_mask, muon_base_mask=muon_base_mask,
          electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
        # to do!
    elif(selectiontype=='chargeflips'):
        # selection of OS events with tight leptons
        # (note: OS/SS selection has to be done in specific selection functions!)
        mask = (ak.sum(ak.concatenate((electron_tight_mask,muon_tight_mask),axis=1),axis=1)==nleptons)
        if(events.metadata['dtype']=='sim'):
            mask = (mask & events.run<0)
    if mask is None:
        msg = 'ERROR in pass_tight_lepton_selection:'
        msg += ' selection type {} not recognized.'.format(selectiontype)
        raise Exception(msg)
    return mask
