#############################
# Test muon reco reweighter #
#############################

# imports
import sys
import os
import argparse
import json
import numpy as np
import ROOT
import array
sys.path.append(os.path.abspath('../../..'))
from plotting.hist2dplotter import plot2dhistogram
import tools.histtools as ht

# input arguments:
parser = argparse.ArgumentParser(description='Plot muon reco scale factors')
parser.add_argument('-s', '--sffile', required=True, type=os.path.abspath)
parser.add_argument('-y', '--year', required=True)
args = parser.parse_args()

# print arguments
print('Running with following configuration:')
for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))

# read json file
with open(args.sffile, 'r') as f:
    sfdata = json.load(f)

# get relevant part of the data
sfdata = sfdata["corrections"][0]
assert sfdata["name"]=="NUM_TrackerMuons_DEN_genTracks"
sfdata = sfdata["data"]

# get eta and pt bins
eta_edges = np.array(sfdata["edges"])
pt_edges = np.array(sfdata["content"][0]["edges"])
shape = (len(eta_edges)-1, len(pt_edges)-1)

# initialize structure
npdata = {
  'nominal': np.zeros(shape),
  'syst': np.zeros(shape),
  'stat': np.zeros(shape) 
}

# loop over eta/pt bins
for eta_idx in range(len(eta_edges)-1):
    for pt_idx in range(len(pt_edges)-1):
        thissfdata = sfdata["content"][eta_idx]['content'][pt_idx]
        # set values
        npdata['nominal'][eta_idx,pt_idx] = thissfdata['content'][0]
        npdata['stat'][eta_idx,pt_idx] = thissfdata['content'][1]
        npdata['syst'][eta_idx,pt_idx] = thissfdata['content'][2]
      
# convert to TH2
nominal = npdata['nominal']
sys = np.sqrt( np.power(npdata['stat'],2) + np.power(npdata['syst'],2) )
hist = ROOT.TH2D('hist', 'Muon reco scale factors',
                len(eta_edges)-1, array.array('f', list(eta_edges)),
                len(pt_edges)-1, array.array('f', list(pt_edges)))
for i in range(1,len(eta_edges)):
    for j in range(1,len(pt_edges)):
        hist.SetBinContent(i,j, nominal[i-1,j-1])
        hist.SetBinError(i,j, sys[i-1,j-1])

# make a plot
figname = args.sffile.replace('.json','.png')
title = 'Muon reco scale factor for {}'.format(args.year)
plot2dhistogram( hist, figname,
  histtitle=title,
  xtitle='|#eta|',
  ytitle='p_{T} (GeV)') 
