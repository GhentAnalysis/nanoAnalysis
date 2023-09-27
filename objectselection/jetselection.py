################################
# Definition of jet selections #
################################

def jetselection(jets, selectionid=None):
    ### perform jet selection
    # input arguments:
    # - jets: awkward array of jets
    #   (e.g. from events.Jet)
    # - selectionid: selection identifier
    # returns:
    # a boolean mask for jets
    
    # switch between selections
    if( selectionid is None ): return (jets.pt > 0.)
    elif( selectionid=='dummy' ): return jetid_dummy(jets)
    elif( selectionid=='run2ul_default'): return jetid_run2ul_default(jets)
    
    # raise error if selection parameters are invalid
    msg = 'ERROR in jetselection:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)


### dummy jet ID for testing ###

def jetid_dummy(jets):
    selection = (
        (jets.pt > 30.)
        & (abs(jets.eta) < 2.4)
    )
    return selection

### default jet ID for Run 2 UL analyses ###

def jetid_run2ul_default(jets):
    selection = (
        (jets.pt > 25.)
        & (abs(jets.eta) < 2.4)
        & (jets.isTight)
    )
    return selection
