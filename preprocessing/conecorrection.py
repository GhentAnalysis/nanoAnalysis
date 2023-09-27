##########################
# Lepton cone correction #
##########################
# Note: requires the presence of the variable jetPtRatio,
#       which is not in standard NanoAOD and should be added
#       by the preprocessing.leptonvariables module
#       or in the skimming step.

# imports
import sys
import os
import awkward as ak

def electron_cone_correction(events,
      electron_fo_mask=None, 
      electron_tight_mask=None,
      correctionfactor=1.):
    ### add cone-corrected electron pt
    # make sure cone correction cannot be applied multiple times
    # note: best approach with events.metadata does not work
    #       as events.metadata does not seem to be editable after creation.
    #       instead just check for ptNoCone branch
    #key = 'electrons_cone_corrected'
    #if( key in events.metadata.keys() ):
    #    if events.metadata[key]:
    if hasattr(events.Electron, 'ptNoCone'):
            msg = 'ERROR: electrons are already cone-corrected.'
            raise Exception(msg)
    # do the correction
    mask = (electron_fo_mask) & (~electron_tight_mask)
    pt = events.Electron.pt
    ptratio = events.Electron.jetPtRatio
    ptcone = ak.where(mask, correctionfactor*pt/ptratio, pt)
    events.Electron = ak.with_field(events.Electron, pt, where='ptNoCone')
    events.Electron = ak.with_field(events.Electron, ptcone, where='pt')
    # also update key access
    # note: appears to be needed to fully propagate the addition of this variable.
    events['Electron'] = events.Electron

def muon_cone_correction(events,
      muon_fo_mask=None,
      muon_tight_mask=None,
      correctionfactor=1.):
    ### add cone-corrected muon pt
    # make sure cone correction cannot be applied multiple times
    # note: best approach with events.metadata does not work
    #       as events.metadata does not seem to be editable after creation.
    #       instead just check for ptNoCone branch
    #key = 'muons_cone_corrected'
    #if( key in events.metadata.keys() ):
    #    if events.metadata[key]:
    if hasattr(events.Muon, 'ptNoCone'):
            msg = 'ERROR: muons are already cone-corrected.'
            raise Exception(msg)
    # do the correction
    mask = (muon_fo_mask) & (~muon_tight_mask)
    pt = events.Muon.pt
    ptratio = events.Muon.jetPtRatio
    ptcone = ak.where(mask, correctionfactor*pt/ptratio, pt)
    events.Muon = ak.with_field(events.Muon, pt, where='ptNoCone')
    events.Muon = ak.with_field(events.Muon, ptcone, where='pt')
    # also update key access
    # note: appears to be needed to fully propagate the addition of this variable.
    events['Muon'] = events.Muon

def lepton_cone_correction(events,
      electron_fo_mask=None,
      electron_tight_mask=None,
      electron_correctionfactor=1.,
      muon_fo_mask=None,
      muon_tight_mask=None,
      muon_correctionfactor=1.):
    electron_cone_correction(
      events,
      electron_fo_mask=electron_fo_mask,
      electron_tight_mask=electron_tight_mask,
      correctionfactor=electron_correctionfactor
    )
    muon_cone_correction(
      events,
      muon_fo_mask=muon_fo_mask,
      muon_tight_mask=muon_tight_mask,
      correctionfactor=muon_correctionfactor
    )
