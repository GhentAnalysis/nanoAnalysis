#####################################
# Definition of electron selections #
#####################################

def electronselection(electrons, selectionid=None):
    ### perform electron selection
    # input arguments:
    # - electrons: awkward array of electrons
    #   (e.g. from events.Electron)
    # - selectionid: selection identifier
    # returns:
    # a boolean mask for electrons
    
    # switch between selections
    if( selectionid is None ): return (electrons.pt > 0.)
    elif( selectionid=='dummy_loose' ): return electronid_dummy_loose(electrons)
    elif( selectionid=='dummy_tight'): return electronid_dummy_tight(electrons)
    elif( selectionid=='run2ul_loose' ): return electronid_run2ul_loose(electrons)
    elif( selectionid=='ttwloose_fo' ): return electronid_ttwloose_fo(electrons)
    elif( selectionid=='ttwloose_tight' ): return electronid_ttwloose_tight(electrons)
    elif( selectionid=='topmvav2_loose'): return electronid_topmvav2_loose(electrons)
    
    # raise error if selection parameters are invalid
    msg = 'ERROR in electronselection:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)


### dummy electron ID for testing ###

def electronid_dummy_loose(electrons):
    selection = (
        (electrons.pt > 10.)
        & (abs(electrons.eta) < 2.5)
    )
    return selection

def electronid_dummy_tight(electrons):
    selection = (
        electronid_dummy_loose(electrons)
        & (electrons.pt > 30.)
    )
    return selection

### loose selection for Run2 UL TOP lepton MVA based IDs
# references:
# - https://github.com/LukaLambrecht/ewkino/blob/4f5a9908fe4a4b5899671738b1d73193e6a6a16c/objectSelection/ElectronSelector.cc#L20
# - Kirill's AN on the TOP lepton MVA (AN-2022-016)

def electronid_run2ul_loose(electrons):
    selection = (
        (electrons.pt > 10.)
        & (abs(electrons.eta) < 2.5)
        & (abs(electrons.dxy) < 0.05)
        & (abs(electrons.dz) < 0.1)
        & (electrons.sip3d < 8.)
        & (electrons.lostHits < 2)
        & (electrons.miniPFRelIso_all < 0.4)
        & ~( (abs(electrons.eta + electrons.deltaEtaSC)>1.4442)
             & (abs(electrons.eta + electrons.deltaEtaSC)<1.566) )
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

def electronid_ttwloose_fo(electrons):
    selection = (
        electronid_run2ul_loose(electrons)
        & (electrons.convVeto)
        & (electrons.tightCharge==2)
        & ( (electrons.mvaTOP > 0.81)
            | ( (electrons.jetBTagDeepFlavor < 0.1)
                & (electrons.jetPtRatio > 0.4)
                & electrons.mvaFall17V2noIso_WPL ) )
    )
    return selection

def electronid_ttwloose_tight(electrons):
    selection = (
        electronid_ttwloose_fo(electrons)
        & (electrons.mvaTOP > 0.81)
    )
    return selection

### loose selection for TOPMVA-v2 preselection
# references:
# - Kirill's AN on the TOP lepton MVA (AN-2022-016)

def electronid_topmvav2_loose(electrons):
    selection = (
        (electrons.isPFcand)
        & (electrons.pt > 10.)
        & (abs(electrons.eta) < 2.5)
        & (electrons.sip3d < 15.)
        & (electrons.miniPFRelIso_all < 1.)
        & ~( (abs(electrons.eta + electrons.deltaEtaSC)>1.4442)
             & (abs(electrons.eta + electrons.deltaEtaSC)<1.566) )
    )
    return selection
