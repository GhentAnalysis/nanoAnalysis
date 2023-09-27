#####################################################
# Functions for cleaning objects around one another #
#####################################################

# imports
import sys
from pathlib import Path
import awkward as ak


def clean(toclean, cleanfrom, conesize):
    ### internal helper function
    # references for use of 'nearest' function:
    # - https://coffeateam.github.io/coffea/notebooks/nanoevents.html
    # - https://github.com/CoffeaTeam/coffea/blob/
    #   58ad8d138905e186e0e0a57239f77026f9a1229e/coffea/nanoevents/methods/vector.py#L706
    dr = toclean.nearest(cleanfrom, return_metric=True)[1]
    mask = (dr > conesize)
    mask = ak.fill_none(mask, True)
    return mask

def clean_electrons_from_muons(electrons, muons, conesize = 0.05):
    ### define a mask to ignore electrons in a small cone around loose muons
    # input arguments:
    # - electrons: awkward array of electrons (e.g. from events.Electron)
    # - muons: awkward array of muons (e.g. from events.Muon)
    #   note: array of muons must match array of electrons at event axis
    #         (e.g. by coming from the same events array)
    #   note: to select only specific muons to clean from,
    #         use e.g. events.Muon[some_mask] instead of events.Muon
    # - conesize: size of the cone around the muon in which to ignore electrons
    # returns: a mask for electrons that do not overlap with selected muons
    return clean(electrons, muons, conesize)

def clean_jets_from_leptons(jets, leptons, conesize=0.4):
    ### define a mask to ignore jets in a cone around leptons
    # input arguments:
    # - jets: awkward array of jets (e.g. from events.Jet)
    # - leptons: awkward array of leptons
    #   note: array of leptons must match array of jets at event axis
    #         (e.g. by coming from the same events array)
    #   note: array of leptons can be constructed from electrons and muons
    #         with e.g. ak.with_name(ak.concatenate((el,mu),axis=1), 'PtEtaPhiMCandidate')
    #   note: to select only specific leptons to clean from,
    #         both electron and muon arrays can be masked before concatenation
    # - conesize: size of the cone around the lepton in which to ignore jets
    # returns: a mask for jets that do not overlap with selected leptons
    return clean(jets, leptons, conesize)
