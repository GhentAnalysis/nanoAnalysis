##########################################################################
# Add custom variables to lepton objects related to generator-level info #
##########################################################################

# imports
import sys
import os
import awkward as ak

def geometricmatch(events, recopart, recopartpdgid, allowphoton=False):
    ### internal helper function to determine geometric gen match,
    # in case the builtin gen matching is not valid.
    # based on:
    # https://github.com/LukaLambrecht/ewkino/blob/nanoaod/objects/src/LeptonGeneratorInfo.cc
    # which was in turn based on:
    # https://github.com/GhentAnalysis/heavyNeutrino/blob/UL_master/multilep/src/GenTools.cc
    # input arguments:
    # - events: object of type NanoEventsArray
    # - recopart: array of reco particles, obtained via e.g. events.Electron
    # - recopartpdgid: pdg id (in absolute value) of reco particles
    # - allowphoton: whether to allow a matching gen particle with pdgid of a photon
    # - verbose: whether to do some info printing
    # returns:
    # array of matching gen particles (with None elements where no valid match was found)
    
    # define mask for which gen particles to consider for matching
    statusmask = (events.GenPart.status==1) # (for stable gen particles)
    if abs(recopartpdgid)==15:
        statusmask = ((events.GenPart.status==2)
                      & events.GenPart.hasFlags(['isLastCopy']))
                      # (special case for taus)
    pdgidmask = (abs(events.GenPart.pdgId)==abs(recopartpdgid))
    if allowphoton: pdgidmask = pdgidmask | (abs(events.GenPart.pdgId)==22)
    genmask = statusmask & pdgidmask

    # find nearest selected gen particle for each reco particle
    (nearestgen,dr) = recopart.nearest(events.GenPart[genmask], threshold=0.2, return_metric=True)

    # repeat the procedure with allowing photons for cases without valid match so far
    if not allowphoton:
        nearestgen2 = geometricmatch(events, recopart, recopartpdgid, allowphoton=True)
        nearestgen = ak.where(ak.is_none(nearestgen, axis=1), nearestgen2, nearestgen)
    
    # return the result
    return nearestgen

def findmatch(events, recopart, recopartpdgid):
    ### internal helper function to determine gen match.
    # priority is given on builtin matching, with fallback to geometric matching.
    # input arguments:
    # - events: object of type NanoEventsArray
    # - recopart: array of reco particles, obtained via e.g. events.Electron
    # - recopartpdgid: pdg id (in absolute value) of reco particles
    # returns:
    # array of matching gen particles (with None elements where no valid match was found)
    genmatch = recopart.matched_gen # (builtin matching)
    validmatch = (genmatch.pdgId==recopart.pdgId)
    validmatch = ak.fill_none(validmatch, False)
    genmatch_geometric = geometricmatch(events, recopart, recopartpdgid)
    genmatch = ak.where(validmatch, genmatch, genmatch_geometric)
    # implementation note: the above is inefficient as geometric matches are calculated
    # for all reco particles, while we only need the ones for which the builtin match is invalid.
    # so ideally would do something like this:
    # genmatch[~validmatch] = geometricmatch(events, recopart[~validmatch], recopartpdgid)
    # however, this kind of in-place assignment seems to be prohibited in awkward arrays.
    return genmatch

def genpart_is_prompt(genpart):
    ### internal helper function to determine promptness of gen particles
    # input arguments:
    # - genpart: array of gen particles, obtained via e.g. events.Electron.matched_gen
    # returns:
    # mask for promptness of the provided gen particles
    isprompt = ((genpart.hasFlags(['isPrompt']))
                | (genpart.hasFlags(['isDirectPromptTauDecayProduct']))
                | (genpart.hasFlags(['isHardProcess']))
                | (genpart.hasFlags(['fromHardProcess']))
                | (genpart.hasFlags(['fromHardProcessBeforeFSR']))
    )
    isprompt = ak.fill_none(isprompt, False)
    return isprompt

def add_electron_gen_variables(events, variables=['all']):
    ### add electron generator variables to events
    # input arguments:
    # - events: an object of type NanoEventsArray
    # - variables: list of variable names to add
    #   (default: all variables defined in this function)
    allvariables = ('all' in variables)
    checkvars = variables[:]
    # do matching one time only (instead of repeating for all variables)
    matches = findmatch(events, events.Electron, 11)
    # isPrompt
    if( allvariables or 'isPrompt' in variables):
        if not allvariables: checkvars.remove('isPrompt')
        events.Electron = ak.with_field(events.Electron, 
          genpart_is_prompt(matches),
          where='isPrompt')
    # matchPdgId
    if( allvariables or 'matchPdgId' in variables):
        if not allvariables: checkvars.remove('matchPdgId')
        events.Electron = ak.with_field(events.Electron,
          ak.fill_none(matches.pdgId, 0),
          where='matchPdgId')
    # isChargeFlip
    if( allvariables or 'isChargeFlip' in variables ):
        if not allvariables: checkvars.remove('isChargeFlip')
        events.Electron = ak.with_field(events.Electron,
          ak.fill_none(events.Electron.pdgId==-matches.pdgId, False),
          where='isChargeFlip')
    # also update key access
    events['Electron'] = events.Electron
    # check unrecognized variables
    if( not allvariables and len(checkvars)>0 ):
        msg = 'ERROR in leptongenvariables.py / add_electron_gen_variables:'
        msg += ' following variables were not recognized: {}.'.format(checkvars)
        raise Exception(msg)

def add_muon_gen_variables(events, variables=['all']):
    ### add muon variables to events
    # input arguments:
    # - events: an object of type NanoEventsArray
    # - variables: list of variable names to add
    #   (default: all variables defined in this function)
    allvariables = ('all' in variables)
    checkvars = variables[:]
    # do matching one time only (instead of repeating for all variables)
    matches = findmatch(events, events.Muon, 13)
    # isPrompt
    if( allvariables or 'isPrompt' in variables):
        if not allvariables: checkvars.remove('isPrompt')
        events.Muon = ak.with_field(events.Muon,
          genpart_is_prompt(matches),
          where='isPrompt')
    # matchPdgId
    if( allvariables or 'matchPdgId' in variables):
        if not allvariables: checkvars.remove('matchPdgId')
        events.Muon = ak.with_field(events.Muon,
          ak.fill_none(matches.pdgId, 0),
          where='matchPdgId')
    # isChargeFlip
    if( allvariables or 'isChargeFlip' in variables ):
        if not allvariables: checkvars.remove('isChargeFlip')
        events.Muon = ak.with_field(events.Muon,
          ak.fill_none(events.Muon.pdgId==-matches.pdgId, False),
          where='isChargeFlip')
    # also update key access
    events['Muon'] = events.Muon
    # check unrecognized variables
    if( not allvariables and len(checkvars)>0 ):
        msg = 'ERROR in leptongenvariables.py / add_muon_gen_variables:'
        msg += ' following variables were not recognized: {}.'.format(checkvars)
        raise Exception(msg)
