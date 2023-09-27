########################################################
# Tools for obtaining modified jet collections and MET #
########################################################
# Note: These tools do not work on 'default' NanoAOD
#       as it requires the additional branches added by the
#       jetMetHelperRun2 module from NanoAOD-Tools.
#       See more info here:
#       https://twiki.cern.ch/twiki/bin/viewauth/CMS/NanoAODTools#JME_jetmet_HelperRun2

# Note: It is not very clear which branches are to be used for the nominal case.
#       It can be observed that e.g. Jet.pt is not centered between Jet.pt_jesTotalUp
#       and Jet.pt_jesTotalDown (often both variations are larger than the "nominal").
#       instead, Jet.pt_nom seems to be a good nominal value in between Up and Down.
#       On the other hand, this branch is not available in 'default' NanoAOD,
#       so for consistency the current approach is to use Jet.pt for nominal and
#       e.g. Jet.pt_jesTotalUp * (Jet.pt / Jet.pt_nom) for variations.
#       To be checked if this is ok.

# Note: Same remark as above for MET: instead of MET.pt, MET.T1Smear_pt
#       seems to be the required central value between variations,
#       but scale it back to MET.pt for consistency with "default" NanoAOD samples. 


import awkward as ak


def get_varied_jets(events, variation):
    ### get varied jet collection
    # basically equivalent to events.Jet,
    # but modify the pt and mass branch according to the variation.
    # note: the variation argument is supposed to be of the form
    #       "jes[source][Up or Down]" or "jer[source][Up or Down]"
    # note: the variation "unclustEn[Up or Down]" is also allowed,
    #       but this will return the nominal jet collection
    jets = events.Jet
    if variation=='nominal': return jets
    if variation.startswith('unclustEn'): return jets
    pt_ratio = jets.pt / jets.pt_nom
    mass_ratio = jets.mass / jets.mass_nom
    jets_pt = getattr(jets, 'pt_{}'.format(variation)) * pt_ratio
    jets_mass = getattr(jets, 'mass_{}'.format(variation)) * mass_ratio
    jets = ak.with_field(jets, jets_pt, where='pt')
    jets = ak.with_field(jets, jets_mass, where='mass')
    return jets

def get_varied_met(events, variation):
    ### equivalent to get_varied_jets but for MET
    # note: the variation argument is supposed to be of the form
    #       "jes[source][Up or Down]" or "jer[source][Up or Down]"
    #       or "unclustEn[Up or Down]"
    met = events.MET
    if variation=='nominal': return met
    pt_ratio = met.pt / met.T1Smear_pt
    phi_diff = met.phi - met.T1Smear_phi
    met_pt = getattr(met, 'T1Smear_pt_{}'.format(variation)) * pt_ratio
    met_phi = getattr(met, 'T1Smear_phi_{}'.format(variation)) + phi_diff
    met = ak.with_field(met, met_pt, where='pt')
    met = ak.with_field(met, met_phi, where='phi')
    return met

def get_available_jec_variations(events):
    ### get list of available JEC variations 
    variations = ([b[3:-2] for b in events.Jet.fields 
                   if (b.startswith('pt_jes') and b.endswith('Up'))])
    return variations

def get_available_jer_variations(events):
    ### get list of available JER variations 
    variations = ([b[3:-2] for b in events.Jet.fields
                   if (b.startswith('pt_jer') and b.endswith('Up'))])
    return variations
