#############################################################
# Tools for working with uncertainties in combine datacards #
#############################################################

# import python modules
import sys
import os
import fnmatch
import math


def get_lumi_unc_dict():
  # note: values to update with latest recommendations!
  res = {}
  res['2016PreVFP'] = {'lumi':1.006,
		       'lumi_2016':1.01}
  res['2016PostVFP'] = res['2016PreVFP']
  res['2017'] = {'lumi':1.009,
		 'lumi_2017':1.02,
		 'lumi_20172018':1.006}
  res['2018'] = {'lumi':1.02,
		 'lumi_2018':1.015,
		 'lumi_20172018':1.002}
  # add root sum square 
  # (only approximate, for naive run2 plotting)
  rss = 0
  for sources in res.values():
    for val in sources.values():
      rss += math.pow(val-1,2)
  rss = 1+math.sqrt(rss)
  res['run2'] = {'lumi': rss}
  return res


def get_lumi_unc_sources( year ):
  return get_lumi_unc_dict()[year]


def get_systematics_to_disable( processes, pnonorm=None, year=None, allyears=None ):
  ### get systematics to disable in a ProcessInfoCollection
  
  rmforall = []
  rmspecific = {}
  for p in processes: rmspecific[p] = []

  # remove overlap between renormalization and factorization scales
  rmforall.append('rScale*')
  rmforall.append('fScale*')
  rmforall.append('rfScales*')
  #rmforall.append('qcdScales*')

  # remove overlap between pdf envelope types
  rmforall.append('pdfShapeEnv*')
  #rmforall.append('pdfShapeRMS*')

  # remove norm uncertainties for specified processes
  if pnonorm is not None:
    for p in pnonorm:
      if p in processes: 
        rmspecific[p].append('rScaleNorm*')
        rmspecific[p].append('fScaleNorm*')
        rmspecific[p].append('rfScalesNorm*')
        rmspecific[p].append('qcdScalesNorm*')
        rmspecific[p].append('pdfNorm*')
        rmspecific[p].append('isrNorm*')
        rmspecific[p].append('fsrNorm*')

  # remove nJets/nBJets uncertainties for all but WZ and ZZ
  for p in processes:
    if( p=='WZ' or p=='ZZ' ): continue
    rmspecific[p].append('njets')
    rmspecific[p].append('nbjets')

  # remove specific nJets uncertainty except for chargeflips
  # (also remove for chargeflips since it was not yet correctly initialized)
  for p in processes:
    for y in allyears: rmspecific[p].append('njetscf{}'.format(y))

  # remove individual qcd and pdf variations
  # (if not done so before)
  rmforall.append('qcdScalesShapeVar*')
  rmforall.append('pdfShapeVar*')

  # remove overlap between JEC sources
  #rmforall.append('JEC')
  rmforall.append('JECGrouped*')
  rmforall.append('JECGrouped_Total*')

  # remove grouped JEC sources for nonprompt
  # (they are not yet correctly initialized)
  # (now commented out since will use single JEC for now)
  #for p in ['Nonprompt', 'NonpromptMu', 'NonpromptE']:
  #  if p in processes: rmspecific[p].append('JECGrouped*')

  # remove fsrShape for WZ
  # (gives unresolved strange behaviour in latest iteration)
  for p in processes:
    if p=='WZ': rmspecific[p].append('fsrShape')

  # remove uncertainties for other years
  # (cleaner, but not strictly needed since they are set to nominal anyway)
  if( year is not None and allyears is not None and year in allyears ):
    for y in allyears:
      if y==year: continue
      rmforall.append('*{}'.format(y))

  return (rmforall, rmspecific)

def remove_systematics_default( processinfo, year=None ):
  ### default sequence removing some shape systematics

  # define processes for which normalization systematics should be removed
  # (because they will get dedicated normalization uncertainties later on)
  pnonorm = ['WZ','ZZ','TTZ','ZG']
  if 'Nonprompt' in processinfo.plist: pnonorm.append('Nonprompt')
  if 'NonpromptE' in processinfo.plist: pnonorm.append('NonpromptE')
  if 'NonpromptMu' in processinfo.plist: pnonorm.append('NonpromptMu')
  if 'Chargeflips' in processinfo.plist: pnonorm.append('Chargeflips')
  for p in processinfo.plist:
    for tag in ['TTW','TTW0','TTW1','TTW2','TTW3','TTW4']:
      if p.startswith(tag): pnonorm.append(p)
  allyears = ['2016PreVFP', '2016PostVFP', '2017', '2018']
  # get standard systematics to disable
  (rmforall, rmspecific) = get_systematics_to_disable( processinfo.plist, 
    pnonorm=pnonorm, year=year, allyears=allyears )
  removedforall = []
  removedspecific = {}
  # remove some systematics for all processes
  for pattern in rmforall:
    for s in processinfo.slist:
      if fnmatch.fnmatch(s, pattern):
        processinfo.disablesys( s, processinfo.plist )
        removedforall.append(s)
  # remove some systematics for selected processes
  for p,patterns in rmspecific.items():
    removedspecific[p] = []
    for pattern in patterns:
      for s in processinfo.slist:
        if fnmatch.fnmatch(s, pattern):
          processinfo.disablesys( s, [p] )
          removedspecific[p].append(s)
  return (removedforall, removedspecific)

def remove_systematics_all( processinfo ):
  ### disable all shape systematics in a processinfo
  for s in processinfo.slist:
    processinfo.disablesys( s, processinfo.plist )
  return (processinfo.slist, dict())

def add_systematics_default( processinfo, year=None ):
  ### default sequence of adding some norm systematics

  # initialize list of all norm systematics
  normsyslist = []

  # add luminosity uncertainty
  if year is not None:
    for source,impact in get_lumi_unc_sources(year).items():
      impacts = {}
      # set impact for all processes
      for p in processinfo.plist: impacts[p] = impact
      processinfo.addnormsys( source, impacts )
      # disable for data-driven backgrounds
      to_disable = ([
        'Nonprompt', 
        'NonpromptE', 'NonpromptMu',
        'Chargeflips'
      ])
      for p in to_disable:
          if p in processinfo.plist:
              processinfo.disablesys( source, [p] )
      normsyslist.append(source)

  # add trigger uncertainty
  # now implemented as reweighter + non-flat uncertainty!
  '''if year is not None:
    impacts = {}
    for p in processinfo.plist: impacts[p] = 1.02
    source = 'Trigger_{}'.format(year)
    processinfo.addnormsys( source, impacts )
    if 'Nonprompt' in processinfo.plist:
      processinfo.disablesys( source, ['Nonprompt'] )
    if 'Chargeflips' in processinfo.plist:
      processinfo.disablesys( source, ['Chargeflips'] )
    normsyslist.append(source)'''

  # add individual norm uncertainties
  norms = ({
    #'WZ': 1.1, now freely floating, rate parameter
    #'ZZ': 1.1, now freely floating, rate parameter
    #'TTZ': 1.15, now freely floating, rate parameter
    'ZG': 1.3, # synced with Oviedo
    'TTG': 1.3, # synced with Oviedo
    'Multiboson': 1.5 # synced with Oviedo
  })
  if 'Nonprompt' in processinfo.plist:
    norms['Nonprompt'] = 1.3
    # (use 30% as equivalent to 20% correlated + 20% uncorrelated in split case)
  if 'NonpromptE' in processinfo.plist:
    norms['NonpromptE'] = 1.2 # contribution for electron nonprompts only
  if 'NonpromptMu' in processinfo.plist:
    norms['NonpromptMu'] = 1.2 # contribution for muon nonprompts only
  if 'Chargeflips' in processinfo.plist:
    norms['Chargeflips'] = 1.3
  for process,mag in norms.items():
    source = 'Norm_{}'.format(process)
    impacts = {}
    for p in processinfo.plist: impacts[p] = '-'
    impacts[process] = mag
    processinfo.addnormsys( source, impacts )
    normsyslist.append(source)

  # also add a correlated uncertainty for split nonprompt
  if( 'NonpromptE' in processinfo.plist and 'NonpromptMu' in processinfo.plist ):
    source = 'Norm_Nonprompt'
    impacts = {}
    for p in processinfo.plist: impacts[p] = '-'
    impacts['NonpromptE'] = 1.2
    impacts['NonpromptMu'] = 1.2
    processinfo.addnormsys( source, impacts )
    normsyslist.append(source)
 
  # return list of all norm systematics that were added
  return normsyslist

def add_systematics_dummy( processinfo ):
  normsyslist = []
  norms = dict()
  for p in processinfo.plist:
    if( p=='Nonprompt' or p=='Chargeflips' ): norms[p] = 1.3
    elif( p.startswith('TTW') ): continue
    else: norms[p] = 1.15
  for process,mag in norms.items():
    source = 'Norm_{}'.format(process)
    impacts = {}
    for p in processinfo.plist: impacts[p] = '-'
    impacts[process] = mag
    processinfo.addnormsys( source, impacts )
    normsyslist.append(source)
  return normsyslist
