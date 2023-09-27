################################################################
# Collection of mappings between code strings and plot strings #
################################################################

def get_region_dict():
  ### get a printable name for each selection region
  res = {}
  res['signalregion_dilepton_inclusive'] = 'SR 2L'
  res['signalregion_dilepton_ee'] = 'SR 2L (ee)'
  res['signalregion_dilepton_em'] = 'SR 2L (em)'
  res['signalregion_dilepton_me'] = 'SR 2L (me)'
  res['signalregion_dilepton_mm'] = 'SR 2L (mm)'
  res['signalregion_dilepton_plus'] = 'SR 2L (++)'
  res['signalregion_dilepton_minus'] = 'SR 2L (--)'
  res['signalregion_trilepton'] = 'SR 3L'
  res['wzcontrolregion'] = 'WZ CR'
  res['zzcontrolregion'] = 'ZZ CR'
  res['zgcontrolregion'] = 'ZG CR'
  res['trileptoncontrolregion'] = '3L CR'
  res['fourleptoncontrolregion'] = '4L CR'
  res['npcontrolregion_met_dilepton_inclusive'] = 'NP CR 2L'
  res['npcontrolregion_met_dilepton_ee'] = 'NP CR 2L (ee)'
  res['npcontrolregion_met_dilepton_em'] = 'NP CR 2L (em)'
  res['npcontrolregion_met_dilepton_me'] = 'NP CR 2L (me)'
  res['npcontrolregion_met_dilepton_mm'] = 'NP CR 2L (mm)'
  res['npcontrolregion_lownjets_dilepton_inclusive'] = 'NP-jets CR 2L'
  res['cfcontrolregion_inclusivejets'] = 'CF CR'
  res['cfcontrolregion_highnjets'] = 'CF+jets CR'
  return res

def get_process_dict():
  ### get a printable name for each process
  res = {}
  res['NonpromptE'] = 'Nonprompt (e)'
  res['NonpromptMu'] = 'Nonprompt (mu)'
  return res
