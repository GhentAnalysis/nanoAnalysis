###################################################################
# Perform TTW event selection and calculate interesting variables #
###################################################################
# input: a single NanoAOD file.
#   it is recommended to not use standard NanoAOD,
#   but rather files preprocessed by this skimmer:
#   https://github.com/GhentAnalysis/nanoSkimming.
# output: a single ROOT file with one or multiple trees.
#   the output ROOT file has the following folder structure:
#   <event selection>/<selection type>/<selection systematic>/Events;
#   the branches of each Events tree hold per-event scalar variables.


# import python modules
import sys
import os
import argparse
from pathlib import Path
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
# import framework modules
sys.path.append(str(Path(__file__).parents[1]))
from objectselection.electronselection import electronselection
from objectselection.muonselection import muonselection
from objectselection.jetselection import jetselection
from objectselection.bjetselection import bjetselection
from objectselection.cleaning import clean_electrons_from_muons
from objectselection.cleaning import clean_jets_from_leptons
from objectselection.conecorrectionfactors import electron_cone_correction_factor
from objectselection.conecorrectionfactors import muon_cone_correction_factor
import objectselection.jetuncs as jetuncs
from preprocessing.conecorrection import lepton_cone_correction
from preprocessing.preprocessor import PreProcessor
from samples.sample import year_from_sample_name
from samples.sample import dtype_from_sample_name
from samples.sampleweights import SampleWeights
import eventselection.lepton_selection_tools as lst
import eventselection.sample_selection_tools as sst
import eventselection.trigger_selection_tools as tst
import tools.argparsetools as apt
from tools.readfakeratetools import readfrmapfromfile
from tools.readchargefliptools import readcfmapfromfile
from reweighting.implementation.run2ulreweighter import get_run2ul_reweighter
# import local modules
sys.path.append(os.path.abspath('eventselections'))
from eventselections import pass_event_selection
sys.path.append(os.path.abspath('eventvariables'))
from eventvariables import calculate_event_variables
sys.path.append(os.path.abspath('systematics'))
from systematics_type import systematics_type
from systematics_tools import get_selection_systematics


if __name__=='__main__':

  sys.stderr.write('###starting###\n')

  # input arguments:
  parser = argparse.ArgumentParser(description='Perform event selection and calculate analysis variables')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
  parser.add_argument('-n', '--nentries', type=int, default=-1)
  parser.add_argument('-s', '--eventselection', required=True, nargs='+')
  parser.add_argument('-t', '--selectiontype', default=['tight'], nargs='+')
  parser.add_argument('--systematics', default=[], nargs='+')
  parser.add_argument('--elfrmap', default=None, type=apt.path_or_none)
  parser.add_argument('--mufrmap', default=None, type=apt.path_or_none)
  parser.add_argument('--elcfmap', default=None, type=apt.path_or_none)
  parser.add_argument('--bdt', default=None, type=apt.path_or_none)
  parser.add_argument('--forcenentries', default=False, action='store_true')
  parser.add_argument('--skimmed', default=False, action='store_true')
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))
  sys.stdout.flush()

  # get sample metadata
  year = year_from_sample_name(args.inputfile)
  dtype = dtype_from_sample_name(args.inputfile)
  samplename = os.path.basename(args.inputfile)
  with uproot.open(args.inputfile) as f:
    events = f['Events']
    nevents = events.num_entries
  print('File paramters:')
  print('  - year {}'.format(year))
  print('  - dtype {}'.format(dtype))
  print('  - available events: {}'.format(nevents))

  # manage systematics
  if dtype=='data': systematics = []
  else:
    if( len(args.systematics)==1 and args.systematics[0]=='all' ):
      args.systematics = list(systematics_type.keys())
      print('Systematics:')
      for systematic in args.systematics: print('  - {}'.format(systematic))
    else:
      for systematic in args.systematics:
        if systematic not in systematics_type.keys():
          raise Exception('ERROR: systematic {} not recognized.'.format(systematic))

  # manage number of entries
  if args.nentries>0:
    if dtype=='data':
      if not args.forcenentries:
        print('WARING: partial file processing for data does not make sense,')
        print('        processing full file instead.')
        args.nentries = -1
      else:
        print('WARING: partial file processing for data does not make sense,')
        print('        use this only for testing.')

  # make NanoEvents array
  print('Loading events from input file...')
  year = year_from_sample_name(args.inputfile)
  dtype = dtype_from_sample_name(args.inputfile)
  samplename = os.path.basename(args.inputfile)
  events = NanoEventsFactory.from_root(
    args.inputfile,
    entry_stop=args.nentries if args.nentries>=0 else None,
    schemaclass=NanoAODSchema,
    metadata={'year': year, 'samplename': samplename, 'dtype': dtype}
  ).events()
  sys.stdout.flush()

  # manage number of entries reweighting
  nentries_reweight = 1
  if args.nentries>0:
    nentries_reweight = nevents / min(args.nentries, nevents)
    print('Using reweighting factor {} because of partial file processing'.format(nentries_reweight))
  nevents = min(args.nentries, nevents)

  # make sample generator weights
  weights = None
  if( dtype=='sim' ): weights = SampleWeights(args.inputfile)

  # load fake rate maps if needed
  electronfrmap = None
  if args.elfrmap is not None:
    electronfrmap = readfrmapfromfile(args.elfrmap, year, 'electron', verbose=True)
  muonfrmap = None
  if args.mufrmap is not None:  
    muonfrmap = readfrmapfromfile(args.mufrmap, year, 'muon', verbose=True)

  # load charge flip maps if needed
  electroncfmap = None
  if args.elcfmap is not None:
    electroncfmap = readcfmapfromfile(args.elcfmap, year, 'electron', verbose=True)

  # initialize output structure
  output_trees = {}
  for eventselection in args.eventselection:
    output_trees[eventselection] = {}
    for selectiontype in args.selectiontype:
      output_trees[eventselection][selectiontype] = {}

  # calculate additional variables
  if args.skimmed: pass
  else:
    preprocessor = PreProcessor()
    leptongenvariables = ['isPrompt'] if dtype=='sim' else []
    preprocessor.process(events,
        leptongenvariables=leptongenvariables,
        leptonvariables=[
          'jetPtRatio',
          'jetBTagDeepFlavor'
        ],
        topmvavariable='mvaTOP', topmvaversion='ULv1',
        dotriggers=True
    )

  # calculate object masks
  print('Performing lepton selection...')
  sys.stdout.flush()
  muon_loose_mask_nominal = muonselection(events.Muon, selectionid='run2ul_loose')
  muon_fo_mask_nominal = muonselection(events.Muon, selectionid='ttwloose_fo')
  muon_tight_mask_nominal = muonselection(events.Muon, selectionid='ttwloose_tight')
  electron_cleaning_mask = clean_electrons_from_muons(events.Electron, events.Muon[muon_loose_mask_nominal])
  electron_loose_mask_nominal = ( 
    (electronselection(events.Electron, selectionid='run2ul_loose'))
    & electron_cleaning_mask )
  electron_fo_mask_nominal = ( 
    (electronselection(events.Electron, selectionid='ttwloose_fo'))
    & electron_cleaning_mask )
  electron_tight_mask_nominal = ( 
    (electronselection(events.Electron, selectionid='ttwloose_tight'))
    & electron_cleaning_mask )
  leptonsforcleaningjets = ak.with_name(ak.concatenate(
    (events.Electron[electron_fo_mask_nominal], events.Muon[muon_fo_mask_nominal]), axis=1),
    'PtEtaPhiMCandidate')
  jet_cleaning_mask = clean_jets_from_leptons(events.Jet, leptonsforcleaningjets)
  jet_mask_nominal = (
    jetselection(events.Jet, selectionid='run2ul_default')
    & jet_cleaning_mask )
  bjet_loose_mask_nominal = (
    jet_mask_nominal 
    & bjetselection(events.Jet, year=year, algo='deepflavor', level='loose') )

  # do lepton cone correction
  lepton_cone_correction(events,
    electron_fo_mask=electron_fo_mask_nominal, electron_tight_mask=electron_tight_mask_nominal,
    muon_fo_mask=muon_fo_mask_nominal, muon_tight_mask=muon_tight_mask_nominal,
    electron_correctionfactor=electron_cone_correction_factor('ttwloose'),
    muon_correctionfactor=muon_cone_correction_factor('ttwloose')
  )

  # make a reweighter
  if dtype=='sim':
    print('Initializing reweighter')
    reweighter = get_run2ul_reweighter(year, dobtagnormalize=True)
    # note: perhaps this should be done inside the loop over selection systematics,
    #       for each selection systematic separately.
    print('Normalizing b-tag reweighter')
    sys.stdout.flush()
    btagreweighter = reweighter.reweighters['btagging']
    jets_for_init = events.Jet[jet_mask_nominal]
    unctypes_for_init = 'all' if 'btagging' in args.systematics else None
    btagreweighter.set_normalization(jets_for_init, unctypes=unctypes_for_init)

  # find systematics for which alternative selections are needed
  if dtype=='data': selection_systematics = {'nominal': ['nominal']}
  else:
    print('Finding systematics that need alternative selections')
    selection_systematics = get_selection_systematics(events, args.systematics, includenominal=True)
    print('Found following selection systematics:')
    for selection_systematic, variations in selection_systematics.items():
      print('  - {}'.format(selection_systematic))
      if( len(variations)>1 or variations[0]!=selection_systematic ):
        for variation in variations: print('    - {}'.format(variation))

  # loop over selection systematics
  electrons_nominal = events.Electron
  muons_nominal = events.Muon
  jets_nominal = events.Jet
  met_nominal = events.MET
  for selection_systematic, variations in selection_systematics.items():
    for variation in variations:
      print('Now running on selection systematic {} ({})'.format(selection_systematic, variation))

      # initialize varied collections and masks to nominal ones
      electrons = electrons_nominal
      muons = muons_nominal
      jets = jets_nominal
      met = met_nominal
      muon_loose_mask = muon_loose_mask_nominal
      muon_fo_mask = muon_fo_mask_nominal
      muon_tight_mask = muon_tight_mask_nominal
      electron_loose_mask = electron_loose_mask_nominal
      electron_fo_mask = electron_fo_mask_nominal
      electron_tight_mask = electron_tight_mask_nominal
      jet_mask = jet_mask_nominal
      bjet_loose_mask = bjet_loose_mask_nominal

      # recalculate some of the objects and masks depending on systematic
      if(selection_systematic=='jec' or selection_systematic=='jer'):
        print('  Recalculating jets and MET')
        jets = jetuncs.get_varied_jets(events, variation)
        events.Jet = jets
        events['Jet'] = jets
        met = jetuncs.get_varied_met(events, variation)
        events.MET = met
        events['MET'] = met
        jet_cleaning_mask_ = clean_jets_from_leptons(events.Jet, leptonsforcleaningjets)
        jet_mask = (
          jetselection(events.Jet, selectionid='run2ul_default')
          & jet_cleaning_mask )
        bjet_loose_mask = (
          jet_mask
          & bjetselection(events.Jet, year=year, algo='deepflavor', level='loose') )
      elif(selection_systematic=='uncl'):
        print('  Recalculating MET')
        met = jetuncs.get_varied_met(events, variation)
        events.MET = met
        events['MET'] = met
      elif(selection_systematic=='nominal'): pass
      else:
        msg = 'ERROR: selection systematic {} not recognized.'.format(selection_systematic)
        raise Exception(msg)
        
      # calculate event variables
      print('  Calculating event variables')
      sys.stdout.flush()
      variables = calculate_event_variables(events,
        weights=weights, nentries_reweight=nentries_reweight,
        electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
        electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask,
        jet_mask=jet_mask, bjet_mask=bjet_loose_mask,
        electronfrmap=electronfrmap, muonfrmap=muonfrmap,
        electroncfmap=electroncfmap)

      # evaluate the reweighter
      if dtype=='sim':
        reweighter_kwargs = ({
          'electron_mask': electron_tight_mask,
          'muon_mask': muon_tight_mask,
          'jet_mask': jet_mask,
          'bjet_mask': bjet_loose_mask
        })
        # calculate nominal event reweighting factors
        print('  Calculate nominal reweighting factors')
        sys.stdout.flush()
        variables['reweight_nominal'] = reweighter.weights(events, **reweighter_kwargs)
        # only for nominal selection: calculate weight systematics
        if selection_systematic=='nominal':
          print('  Calculate systematic reweighting factors')
          sys.stdout.flush()
          weight_systematics = ([systematic for systematic in args.systematics 
            if systematics_type[systematic]=='weight'])
          for systematic in weight_systematics:
            print('    - {}'.format(systematic))
            sys.stdout.flush()
            unctypes = reweighter.get_uncertainties(systematic)
            if unctypes is None: unctypes = [None]
            for unctype in unctypes:
              key = systematic
              if unctype is not None:
                key += '_'+unctype
                print('      - {}'.format(unctype))
                sys.stdout.flush()
              up = reweighter.weightsup(events, systematic, unctype=unctype, **reweighter_kwargs)
              down = reweighter.weightsdown(events, systematic, unctype=unctype, **reweighter_kwargs)
              variables['reweight_{}_up'.format(key)] = up
              variables['reweight_{}_down'.format(key)] = down

      # loop over event selections and selection types
      for eventselection in args.eventselection:
        for selectiontype in args.selectiontype:
          print('  Performing event selection for {} / {}'.format(
            eventselection, selectiontype))
          sys.stdout.flush()
          # do event selection
          events_mask = pass_event_selection(events, eventselection,
            selectiontype=selectiontype,
            electron_fo_mask=electron_fo_mask, muon_fo_mask=muon_fo_mask,
            electron_tight_mask=electron_tight_mask, muon_tight_mask=muon_tight_mask,
            jet_mask=jet_mask, bjet_mask=bjet_loose_mask)
          print('  Event selection: {} / {}'.format(eventselection, selectiontype))
          print('    Selected {} out of {} events.'.format(ak.sum(events_mask), nevents))
          sys.stdout.flush()
          # make output tree dict
          tree = {}
          for key,val in variables.items():
            tree[key] = val[events_mask]
          systematic_tree_name = selection_systematic
          if( len(variations)>1 or variation!=selection_systematic ):
            systematic_tree_name += '_' + variation
          output_trees[eventselection][selectiontype][systematic_tree_name] = tree

      # reset event objects to nominal
      events.Electron = electrons_nominal
      events['Electron'] = electrons_nominal
      events.Muon = muons_nominal
      events['Muon'] = muons_nominal
      events.Jet = jets_nominal
      events['Jet'] = jets_nominal
      events.MET = met_nominal
      events['MET'] = met_nominal

  # write output trees to file
  outputdir = os.path.dirname(args.outputfile)
  if not os.path.exists(outputdir): os.makedirs(outputdir)
  with uproot.recreate(args.outputfile, compression=uproot.LZMA(9)) as f:
    for eventselection in output_trees.keys():
      for selectiontype in output_trees[eventselection].keys():
        for systematic in output_trees[eventselection][selectiontype].keys():
          print('Writing output for {} / {} / {} ...'.format(
            eventselection, selectiontype, systematic))
          sys.stdout.flush()
          fkey = '{}/{}/{}/Events'.format(eventselection,selectiontype,systematic)
          f[fkey] = output_trees[eventselection][selectiontype][systematic]

  sys.stderr.write('###done###\n')
