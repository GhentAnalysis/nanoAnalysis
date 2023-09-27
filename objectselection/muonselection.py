#################################
# Definition of muon selections #
#################################

def muonselection(muons, selectionid=None):
    ### perform muon selection
    # input arguments:
    # - muons: awkward array of muons
    #   (e.g. from events.Muon)
    # - selectionid: selection identifier
    # returns:
    # a boolean mask for muons
    
    # switch between selections
    if( selectionid is None ): return (muons.pt > 0.)
    elif( selectionid=='dummy_loose' ): return muonid_dummy_loose(muons)
    elif( selectionid=='dummy_tight'): return muonid_dummy_tight(muons)
    elif( selectionid=='run2ul_loose' ): return muonid_run2ul_loose(muons)
    elif( selectionid=='ttwloose_fo' ): return muonid_ttwloose_fo(muons)
    elif( selectionid=='ttwloose_tight' ): return muonid_ttwloose_tight(muons)
    elif( selectionid=='topmvav2_loose'): return muonid_topmvav2_loose(muons)

    # raise error if selection parameters are invalid
    msg = 'ERROR in muonselection:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)


### dummy muon ID for testing ###

def muonid_dummy_loose(muons):
    selection = (
        (muons.pt > 10.)
        & (abs(muons.eta) < 2.4)
    )
    return selection

def muonid_dummy_tight(muons):
    selection = (
        muonid_dummy_loose(muons)
        & (muons.pt > 30.)
    )
    return selection

### loose selection for Run2 UL TOP lepton MVA based IDs
# references:
# - https://github.com/LukaLambrecht/ewkino/blob/4f5a9908fe4a4b5899671738b1d73193e6a6a16c/objectSelection/MuonSelector.cc#L16
# - Kirill's AN on the TOP lepton MVA (AN-2022-016)

def muonid_run2ul_loose(muons):
    selection = (
        (muons.isPFcand)
        & ((muons.isTracker) | (muons.isGlobal))
        & (muons.pt > 10.)
        & (abs(muons.eta) < 2.4)
        & (abs(muons.dxy) < 0.05)
        & (abs(muons.dz) < 0.1)
        & (muons.sip3d < 8.)
        & (muons.miniPFRelIso_all < 0.4)
        & (muons.mediumId)
    )
    return selection

### FO and tight selection for TTTT / TTW (loose) analyses
# references:
# - TTTT AN-2021-182
# - TTW AN-2023-077
# notes:
# - requires the mvaTOP variable to be set, which is not in the standard nanoAOD.
# - also the variables jetPtRatio and jetBTagDeepFlavor are not in standard nanoAOD.
# - so far, only the 2018 FO ID has been implemented (for testing),
#   later extend to other years (with slightly modified values)

def muonid_ttwloose_fo(muons):
    selection = (
        muonid_run2ul_loose(muons)
        & ( (muons.mvaTOP > 0.64)
            | ( (muons.jetBTagDeepFlavor < 0.025)
                & (muons.jetPtRatio > 0.45) ) )
    )
    return selection

def muonid_ttwloose_tight(muons):
    selection = (
        muonid_ttwloose_fo(muons)
        & (muons.mvaTOP > 0.64)
    )
    return selection

### loose selection for TOPMVA-v2 preselection
# references:
# - Kirill's AN on the TOP lepton MVA (AN-2022-016)

def muonid_topmvav2_loose(muons):
    selection = (
        (muons.isPFcand)
        & ((muons.isTracker) | (muons.isGlobal))
        & (muons.pt > 10.)
        & (abs(muons.eta) < 2.4)
        & (muons.sip3d < 15.)
        & (muons.miniPFRelIso_all < 1.)
    )
    return selection
