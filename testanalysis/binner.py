##################################
# Perform event variable binning #
##################################
# input: a single ROOT file processed by eventloop.py.
#   expected input format: a ROOT file with one or multiple trees,
#   named as <event selection>/<selection type>/<selection systematic>/Events;
#   the branches of each tree should hold per-event scalar variables.
# output: a file containing ROOT histograms with binned event variables.
#   output format: a ROOT file with histograms named as follows:
#   <process tag>_<event selection>_<selection_type>_<variable>_<systematic>
#   where the process tag is omitted if not provided as command line arg.


# import python modules
import sys
import os
import argparse
from pathlib import Path
import numpy as np
import awkward as ak
import uproot
import ROOT
# import framework modules
sys.path.append(str(Path(__file__).parents[1]))
from samples.sample import year_from_sample_name
from samples.sample import dtype_from_sample_name
import tools.argparsetools as apt
from tools.variabletools import read_variables


if __name__=='__main__':

  sys.stderr.write('###starting###\n')

  # input arguments:
  parser = argparse.ArgumentParser(description='Perform event binning')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputfile', required=True, type=os.path.abspath)
  parser.add_argument('-n', '--nentries', type=int, default=-1)
  parser.add_argument('-v', '--variables', required=True, type=os.path.abspath)
  parser.add_argument('-s', '--eventselection', default=['auto'], nargs='+')
  parser.add_argument('-t', '--selectiontype', default=['auto'], nargs='+')
  parser.add_argument('--systematics', default=['auto'], nargs='+')
  parser.add_argument('--xsec', default=1, type=float)
  parser.add_argument('--lumi', default=1, type=float)
  parser.add_argument('--process', default=None)
  parser.add_argument('--split', default=False, action='store_true')
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

  # initializations
  year = year_from_sample_name(args.inputfile)
  dtype = dtype_from_sample_name(args.inputfile)
  samplename = os.path.basename(args.inputfile)

  # read variables
  variables = read_variables( args.variables )

  # open the input file
  with uproot.open(args.inputfile) as inputfile:

    # extract the event selections and selection types and selection systematics
    selectionlevels = {}
    if 'auto' in args.eventselection:
      eventselections = ([el.split(';')[0] for el in inputfile.keys()
        if el.count('/')==0])
    else: eventselections = args.eventselection
    for eventselection in eventselections:
      selectionlevels[eventselection] = {}
      if 'auto' in args.selectiontype:
        selectiontypes = ([el.split('/')[1].split(';')[0] for el in inputfile.keys() 
          if (eventselection in el and el.count('/')==1)])
      else: selectiontypes = args.selectiontype
      for selectiontype in selectiontypes:
        if 'auto' in args.systematics:
          selection_variations = ([el.split('/')[2].split(';')[0] for el in inputfile.keys()
            if (eventselection in el and selectiontype in el and el.count('/')==2)])
        else: raise Exception('Not yet implemented')
        selectionlevels[eventselection][selectiontype] = selection_variations

    # extract the weight systematics from nominal tree
    nominaltree = inputfile[eventselection][selectiontype]['nominal']['Events']
    weight_systematics = ([branch.replace('reweight_','').replace('_nom','') 
                           for branch in nominaltree.keys()
                           if(branch.startswith('reweight_') and branch.endswith('_nom'))])
    # (e.g. weight_systematics = ['muonid'])
    weight_variations = ([branch.replace('reweight_','') for branch in nominaltree.keys()
                           if( branch.startswith('reweight_') and not branch.endswith('_nom'))])
    # (e.g. weight_variations = ['nominal', 'muonid_stat_up', 'muonid_stat_down', etc.])
    weight_vartonom = {}
    for weight_variation in weight_variations:
        if weight_variation=='nominal': continue
        match = None
        for weight_systematic in weight_systematics:
            if '{}'.format(weight_systematic) in weight_variation:
                match = weight_systematic
                break
        weight_vartonom[weight_variation] = match

    # printouts of the above
    print('Event selection levels:')
    for eventselection in selectionlevels.keys():
      print('  - {}'.format(eventselection))
      for selectiontype in selectionlevels[eventselection].keys():
        print('    - {}'.format(selectiontype))
        for selection_variation in selectionlevels[eventselection][selectiontype]:
          print('      - {}'.format(selection_variation))
    print('Weight variations:')
    for weight_variation in weight_variations:
      txt = '  - {}'.format(weight_variation)
      if weight_variation!='nominal':
          txt += ' (variation on {})'.format(weight_vartonom[weight_variation])
      print(txt)

    # initialize dict of histograms
    hists = {}

    # loop over event selections and selection types
    for eventselection in selectionlevels.keys():
      hists[eventselection] = {}
      for selectiontype in selectionlevels[eventselection]:
        hists[eventselection][selectiontype] = {}
        for selection_variation in selectionlevels[eventselection][selectiontype]:
          variations = [selection_variation]
          if selection_variation=='nominal':
            if dtype=='sim': variations = weight_variations
            else: variations = ['nominal']
          for variation in variations:
            hists[eventselection][selectiontype][variation] = []
            print('Now running on {} / {} / {}'.format(
              eventselection, selectiontype, variation))
            
            # read the tree
            events = inputfile[eventselection][selectiontype][selection_variation]['Events']
            nevents = events.num_entries
            print('Found tree with {} entries.'.format(nevents))
            # convert to group of arrays
            events = events.arrays(library='np')
            
            # calculate correct event weights (before reweighting)
            if dtype=='data': weights = np.ones(nevents)
            else: weights = events['genNormWeight'] * args.xsec * args.lumi
            if( selectiontype=='fakerate' 
              or selectiontype=='efakerate' 
              or selectiontype=='mfakerate' ):
              frweights = events['fakeRateWeight']
              if dtype=='sim': frweights = -frweights
              weights = np.multiply(weights, frweights)
            if selectiontype=='chargeflips':
              cfweights = events['chargeFlipWeight']
              if dtype=='sim': cfweights = np.zeros(nevents)
              weights = np.multiply(weights, cfweights)
            
            # modify event weights with reweighting factor
            # note: the reweighting factors in the nominal tree
            #       are individual reweighting factors,
            #       so some extra arithmetic is needed to calculate the total weight
            if dtype=='sim':
              nomweights = events['reweight_nominal']
              reweight = nomweights # default case: just use nominal weights
              if( selection_variation=='nominal' ):
                if( variation!='nominal' ):
                  nomweights_single = events['reweight_{}_nom'.format(weight_vartonom[variation])]
                  nomweights_withoutsingle = np.divide( nomweights, nomweights_single )
                  reweight = np.multiply(nomweights_withoutsingle, events['reweight_{}'.format(variation)])
              weights = np.multiply(weights, reweight)
            
            # modify process name if needed
            thisprocess = args.process
            if thisprocess is not None:
              if selectiontype=='fakerate': thisprocess = 'Nonprompt'
              if selectiontype=='efakerate': thisprocess = 'NonpromptE'
              if selectiontype=='mfakerate': thisprocess = 'NonpromptMu'
              if selectiontype=='chargeflips': thisprocess = 'Chargeflips'
            
            # loop over variables
            for variable in variables:
              if variable.variable not in events.keys():
                msg = 'WARNING: variable {} not found in input tree, skipping...'.format(variable.variable)
                print(msg)
                continue
              # initialize a histogram
              histname = '{}_{}_{}_{}'.format(eventselection, selectiontype, variable.name, variation)
              if thisprocess is not None: histname = '{}_{}'.format(thisprocess, histname)
              hist = variable.initialize_histogram(histname=histname)
              # clip the values to be inside the histogram range
              nbins = hist.GetNbinsX()
              minvalue = hist.GetBinLowEdge(1) + hist.GetBinWidth(1)/2.
              maxvalue = hist.GetBinLowEdge(nbins) + hist.GetBinWidth(nbins)/2.
              values = np.clip(events[variable.variable], minvalue, maxvalue)
              # fill the histogram
              for val, weight in zip(values, weights):
                hist.Fill(val, weight)
              hists[eventselection][selectiontype][variation].append(hist)

  # write to output file
  # (note: deliberately using PyROOT instead of uproot
  #  because the latter seems to be unreasonably slow)
  # case 1: single output file for all histograms obtained from input file
  if not args.split:
    outputfile = ROOT.TFile.Open(args.outputfile, 'recreate')
    print('Writing histograms to output file...')  
    # write all histograms
    for eventselection in hists.keys():
      for selectiontype in hists[eventselection].keys():
        for variation in hists[eventselection][selectiontype].keys():
          for hist in hists[eventselection][selectiontype][variation]:
            hist.Write()
    outputfile.Close()
  
  # case 2: multiple output files split per eventselection and selectiontype
  if args.split:
    outputdir = os.path.dirname(args.outputfile)
    print('Writing histograms to output file...')
    for eventselection in hists.keys():
      for selectiontype in hists[eventselection].keys():
        thisoutputdir = os.path.join(outputdir, eventselection, selectiontype)
        if not os.path.exists(thisoutputdir): os.makedirs(thisoutputdir)
        outputfile = os.path.join(thisoutputdir, os.path.basename(args.outputfile))
        outputfile = ROOT.TFile.Open(outputfile, 'recreate')
        for systematic in hists[eventselection][selectiontype].keys():
          for hist in hists[eventselection][selectiontype][systematic]:
            hist.Write()
        outputfile.Close()

  sys.stderr.write('###done###\n')
