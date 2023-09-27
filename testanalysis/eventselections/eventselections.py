##################################
# Definition of event selections #
##################################


# imports
import sys
import os
from pathlib import Path
import awkward as ak
sys.path.append(str(Path(__file__).parents[2]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.jetselection import jetselection
from objectselection.bjetselection import bjetselection
from objectselection.cleaning import clean_electrons_from_muons
from objectselection.cleaning import clean_jets_from_leptons
import eventselection.lepton_selection_tools as lst
import eventselection.sample_selection_tools as sst
import eventselection.trigger_selection_tools as tst
from eventreconstruction.zreco import ZReco
from constants.particlemasses import m_Z


def pass_event_selection(events, eventselection, **kwargs):
    if(eventselection=='signalregion_dilepton_inclusive'):
        return pass_signalregion_dilepton_inclusive(events, **kwargs)
    elif(eventselection=='signalregion_trilepton'):
        return pass_signalregion_trilepton(events, **kwargs)
    elif(eventselection=='trileptoncontrolregion'):
        return pass_trileptoncontrolregion(events, **kwargs)
    elif(eventselection=='fourleptoncontrolregion'):
        return pass_fourleptoncontrolregion(events, **kwargs)
    elif(eventselection=='npcontrolregion_met_dilepton_inclusive'):
        return pass_npcontrolregion_met_dilepton_inclusive(events, **kwargs)
    elif(eventselection=='npcontrolregion_lownjets_dilepton_inclusive'):
        return pass_npcontrolregion_lownjets_dilepton_inclusive(events, **kwargs)
    elif(eventselection=='cfcontrolregion_inclusivejets'):
        return pass_cfcontrolregion_inclusivejets(events, **kwargs)
    elif(eventselection=='cfcontrolregion_highnjets'):
        return pass_cfcontrolregion_highnjets(events, **kwargs)
    else:
        msg = 'ERROR in pass_event_selection:'
        msg += ' event selection {} not recognized.'.format(eventselection)
        raise Exception(msg)

def pass_signalregion_dilepton_inclusive(events,
    selectiontype='tight',
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    cutflow=False):
    # define masks
    metfilter_mask = tst.pass_met_filters(events)
    trigger_mask = tst.pass_any_lepton_trigger(events)
    fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==2)
    lowmass_mask = lst.pass_mll_lowmass_veto(events,
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    photon_mask = sst.pass_photon_overlap_removal(events, samplename=events.metadata['samplename'])
    tight_mask = lst.pass_tight_lepton_selection(events, 2, selectiontype,
      electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
      electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
    # (note: after this point, use electron_fo_mask and muon_fo_mask;
    #  they are equivalent to the tight masks except in case of 'fakerate' selection,
    #  where FO rather than tight leptons are needed in the following steps.)
    pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.),
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    invmass = lst.get_lepton_invmass(events,
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    invmass_mask = (invmass > 30.)
    ss_mask = (abs(ak.sum(ak.concatenate((events.Electron[electron_fo_mask].charge,
      events.Muon[muon_fo_mask].charge),axis=1),axis=1))==2)
    if selectiontype=='chargeflips': ss_mask = ~ss_mask
    zveto_mask = (
      (ak.sum(electron_fo_mask,axis=1)!=2)
      | (abs(invmass - m_Z) > 10.) )
    met_mask = ( events.MET.pt > 30. )
    nbjet_mask = (ak.sum(bjet_mask,axis=1) >= 2)
    njet_mask = (ak.sum(jet_mask,axis=1) >= 3)
    # aggregate masks
    masks = {
      'MET filters': metfilter_mask,
      'Trigger': trigger_mask,
      '2 FO leptons': fo_mask,
      'Low mass veto': lowmass_mask,
      'Photon overlap': photon_mask,
      '2 tight leptons': tight_mask,
      'pT thresholds': pt_mask,
      'Invariant mass veto': invmass_mask,
      'Same sign': ss_mask,
      'Electron Z veto': zveto_mask,
      'MET': met_mask,
      'b-tagged jets': nbjet_mask,
      'Jets': njet_mask
    }
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_signalregion_dilepton_ee(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
    # define masks
    inclusive_masks = pass_signalregion_dilepton_inclusive(
      events, electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
      cutflow=True, **kwargs)
    flavour_mask = (ak.num(events.Electron[electron_fo_mask])==2)
    # aggregate masks
    masks = inclusive_masks
    masks['Lepton flavour (ee)'] = flavour_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_signalregion_dilepton_em(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
    raise Exception('Not yet implemented.')

def pass_signalregion_dilepton_me(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
   raise Exception('Not yet implemented.')

def pass_signalregion_dilepton_mm(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
    # define masks
    inclusive_masks = pass_signalregion_dilepton_inclusive(
      events, electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
      cutflow=True, **kwargs)
    flavour_mask = (ak.num(events.Muon[muon_fo_mask])==2)
    # aggregate masks
    masks = inclusive_masks
    masks['Lepton flavour (mm)'] = flavour_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_signalregion_dilepton_plus(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
    # define masks
    inclusive_masks = pass_signalregion_dilepton_inclusive(
      events, electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
      cutflow=True, **kwargs)
    sign_mask = (ak.sum(ak.concatenate((events.Electron[electron_fo_mask].charge,
      events.Muon[muon_fo_mask].charge),axis=1),axis=1)==2)
    # aggregate masks
    masks = inclusive_masks
    masks['Lepton sign (++)'] = sign_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_signalregion_dilepton_minus(events,
    electron_fo_mask=None, muon_fo_mask=None,
    cutflow=False,
    **kwargs):
    # define masks
    inclusive_masks = pass_signalregion_dilepton_inclusive(
      events, electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
      cutflow=True, **kwargs)
    sign_mask = (ak.sum(ak.concatenate((events.Electron[electron_fo_mask].charge,
      events.Muon[muon_fo_mask].charge),axis=1),axis=1)==-2)
    # aggregate masks
    masks = inclusive_masks
    masks['Lepton sign (--)'] = sign_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_signalregion_trilepton(events,
    selectiontype='tight',
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    cutflow=False):
    # define masks
    metfilter_mask = tst.pass_met_filters(events)
    trigger_mask = tst.pass_any_lepton_trigger(events)
    fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==3)
    lowmass_mask = lst.pass_mll_lowmass_veto(events,
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    photon_mask = sst.pass_photon_overlap_removal(events, samplename=events.metadata['samplename'])
    tight_mask = lst.pass_tight_lepton_selection(events, 3, selectiontype,
      electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
      electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
    pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.,15.),
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    zreco = ZReco(events, halfwindow=10., electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    zveto_mask = ~zreco.has_ztoll_candidate()
    nbjet_mask = (ak.sum(bjet_mask,axis=1) >= 2)
    njet_mask = (ak.sum(jet_mask,axis=1) >= 3)
    # aggregate masks
    masks = {
      'MET filters': metfilter_mask,
      'Trigger': trigger_mask,
      '3 FO leptons': fo_mask,
      'Low mass veto': lowmass_mask,
      'Photon overlap': photon_mask,
      '3 tight leptons': tight_mask,
      'pT thresholds': pt_mask,
      'Z veto': zveto_mask,
      'b-tagged jets': nbjet_mask,
      'Jets': njet_mask
    }
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_trileptoncontrolregion(events,
    selectiontype='tight',
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    cutflow=False):
    # define masks
    metfilter_mask = tst.pass_met_filters(events)
    trigger_mask = tst.pass_any_lepton_trigger(events)
    fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==3)
    lowmass_mask = lst.pass_mll_lowmass_veto(events,
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    photon_mask = sst.pass_photon_overlap_removal(events, samplename=events.metadata['samplename'])
    tight_mask = lst.pass_tight_lepton_selection(events, 3, selectiontype,
      electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
      electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
    pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.,15.),
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    zreco = ZReco(events, halfwindow=10., electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    z_mask = zreco.has_ztoll_candidate()
    # aggregate masks
    masks = {
      'MET filters': metfilter_mask,
      'Trigger': trigger_mask,
      '2 FO leptons': fo_mask,
      'Low mass veto': lowmass_mask,
      'Photon overlap': photon_mask,
      '3 tight leptons': tight_mask,
      'pT thresholds': pt_mask,
      '3 candidate': z_mask
    }
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_fourleptoncontrolregion(events,
    selectiontype='tight',
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    cutflow=False):
    # define masks
    metfilter_mask = tst.pass_met_filters(events)
    trigger_mask = tst.pass_any_lepton_trigger(events)
    fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==4)
    photon_mask = sst.pass_photon_overlap_removal(events, samplename=events.metadata['samplename'])
    tight_mask = lst.pass_tight_lepton_selection(events, 4, selectiontype,
      electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
      electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
    pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.,15.,10.),
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    zreco = ZReco(events, halfwindow=10., electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    z_mask = zreco.has_ztoll_candidate()
    # aggregate masks
    masks = {
      'MET filters': metfilter_mask,
      'Trigger': trigger_mask,
      '4 FO leptons': fo_mask,
      'Photon overlap': photon_mask,
      '4 tight leptons': tight_mask,
      'pT thresholds': pt_mask,
      'Z candidate': z_mask
    }
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_npcontrolregion_met_dilepton_inclusive(events, cutflow=False, **kwargs):
    # define masks
    masks = pass_signalregion_dilepton_inclusive(events, cutflow=True, **kwargs)
    masks['MET'] = ~masks['MET']
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_npcontrolregion_lownjets_dilepton_inclusive(events,
    jet_mask=None, bjet_mask=None, cutflow=False, **kwargs):
    # define masks
    masks = pass_signalregion_dilepton_inclusive(events,
              jet_mask=jet_mask, bjet_mask=bjet_mask, cutflow=True, **kwargs)
    nbjet_mask = (ak.sum(bjet_mask,axis=1) >= 1)
    njet_mask = ((ak.sum(jet_mask,axis=1) >= 1) & (ak.sum(jet_mask,axis=1) <= 2))
    masks['b-tagged jets'] = nbjet_mask
    masks['Jets'] = njet_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_cfcontrolregion_inclusivejets(events,
    selectiontype='tight',
    electron_fo_mask=None, muon_fo_mask=None,
    electron_tight_mask=None, muon_tight_mask=None,
    jet_mask=None, bjet_mask=None,
    cutflow=False):
    # define masks
    metfilter_mask = tst.pass_met_filters(events)
    trigger_mask = tst.pass_any_lepton_trigger(events)
    fo_mask = (ak.sum(ak.concatenate((electron_fo_mask,muon_fo_mask),axis=1),axis=1)==2)
    photon_mask = sst.pass_photon_overlap_removal(events, samplename=events.metadata['samplename'])
    tight_mask = lst.pass_tight_lepton_selection(events, 2, selectiontype,
      electron_base_mask=electron_fo_mask, muon_base_mask=muon_fo_mask,
      electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask)
    pt_mask = lst.pass_lepton_pt_thresholds(events, pt_thresholds=(25.,15.),
      electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    flavour_mask = (ak.sum(electron_fo_mask, axis=1)==2)
    ss_mask = (abs(ak.sum(ak.concatenate((events.Electron[electron_fo_mask].charge,
      events.Muon[muon_fo_mask].charge),axis=1),axis=1))==2)
    if selectiontype=='chargeflips': ss_mask = ~ss_mask
    zreco = ZReco(events, halfwindow=10., samesign=(selectiontype!='chargeflips'),
                  electron_mask=electron_fo_mask, muon_mask=muon_fo_mask)
    z_mask = zreco.has_ztoll_candidate()
    # aggregate masks
    masks = {
      'MET filters': metfilter_mask,
      'Trigger': trigger_mask,
      '2 FO leptons': fo_mask,
      'Photon overlap': photon_mask,
      '2 tight leptons': tight_mask,
      'pT thresholds': pt_mask,
      'Lepton flavour (ee)': flavour_mask,
      'Same sign': ss_mask,
      'Z candidate': z_mask
    }
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask

def pass_cfcontrolregion_highnjets(events,
    jet_mask=None, bjet_mask=None,
    cutflow=False, **kwargs):
    # define masks
    masks = pass_cfcontrolregion_inclusivejets(events,
              jet_mask=jet_mask, bjet_mask=bjet_mask,
              cutflow=True, **kwargs)
    nbjet_mask = (ak.sum(bjet_mask,axis=1) >= 1)
    njet_mask = (ak.sum(jet_mask,axis=1) >= 3)
    masks['b-tagged jets'] = nbjet_mask
    masks['Jets'] = njet_mask
    # return full set of masks for cutflow
    if cutflow: return masks
    # else return only total mask
    maskkeys = list(masks.keys())
    totalmask = masks[maskkeys[0]]
    for key in maskkeys[1:]: totalmask = totalmask & masks[key]
    return totalmask
