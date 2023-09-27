#########################################
# Definition of cone correction factors #
#########################################

def electron_cone_correction_factor(selectionid):
    # switch between selections
    if( selectionid=='dummy' ): return 0.8
    elif( selectionid=='ttwloose' ): return 0.72
    # raise error if selection parameters are invalid
    msg = 'ERROR in electron_cone_correction_factor:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)

def muon_cone_correction_factor(selectionid):
    # switch between selections
    if( selectionid=='dummy' ): return 0.8
    elif( selectionid=='ttwloose' ): return 0.66
    # raise error if selection parameters are invalid
    msg = 'ERROR in muon_cone_correction_factor:'
    msg += ' selection {} not recognized.'.format(selectionid)
    raise Exception(msg)
