###############################
# Definition of b-tagged jets #
###############################

def inbtagacceptance(jets):
    ### internal helper function
    # todo: verify if these selections still correspond
    #       to the maximum acceptance of b-tagging algorithms.
    return ( (jets.pt > 25.) & (abs(jets.eta)<2.4) )

def bjetselection(jets, year=None, algo=None, level=None):
    ### perform b-jet selection
    # input arguments:
    # - jets: awkward array of jets
    #   (e.g. from events.Jet)
    # - year: year identifier for the provided jets
    # - algo: b-tagging algorithm (for now only 'deepflavor' is supported)
    # - level: choose from 'loose', 'medium' or 'tight'
    # returns:
    # a boolean mask for b-tagged jets

    # check arguments
    if( year is None or algo is None or level is None ):
        msg = 'ERROR in bjetselection:'
        msg += ' year, algo and level must all be provided.'
        raise Exception(msg)

    # switch between algorithms
    if( algo=='deepflavor' ): return deepflavorselection(jets, year, level)
    
    # raise error if b-tag algorithm is invalid
    msg = 'ERROR in bjetselection:'
    msg += ' algo {} not recognized.'.format(algo)
    raise Exception(msg)


### DeepJet/DeepFlavor b-tagging ###

def deepflavorselection(jets, year, level):
    ### perform b-jet selection using the DeepJet/DeepFlavor algorithm
    # note: the threshold values can be found here:
    #       https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016preVFP/
    #       (and similar for other years)
    bscores = jets.btagDeepFlavB
    selection = None
    if year=='2016PreVFP':
        if level=='loose': selection = (bscores > 0.0508)
        if level=='medium': selection = (bscores > 0.2598)
        if level=='tight': selection = (bscores > 0.6502)
    if year=='2016PostVFP':
        if level=='loose': selection = (bscores > 0.0480)
        if level=='medium': selection = (bscores > 0.2489)
        if level=='tight': selection = (bscores > 0.6377)
    if year=='2017':
        if level=='loose': selection = (bscores > 0.0532)
        if level=='medium': selection = (bscores > 0.3040)
        if level=='tight': selection = (bscores > 0.7476)
    if year=='2018':
        if level=='loose': selection = (bscores > 0.0490)
        if level=='medium': selection = (bscores > 0.2783)
        if level=='tight': selection = (bscores > 0.7100)
    if selection is None:
        msg = 'ERROR in deepflavorselection:'
        msg += ' year {}, level {} not recognized.'.format(year, level)
        raise Exception(msg)
    acceptance = inbtagacceptance(jets)
    return (acceptance & selection)
