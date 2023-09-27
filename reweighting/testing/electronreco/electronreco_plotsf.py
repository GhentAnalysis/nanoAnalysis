#################################
# Test electron reco reweighter #
#################################

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
parser = argparse.ArgumentParser(description='Plot electron reco scale factors')
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
assert sfdata["name"]=="UL-Electron-ID-SF"
sfdata = sfdata["data"]["content"][0]["value"]["content"]
sfnomdata = sfdata[0]
sfupdata = sfdata[1]
sfdowndata = sfdata[2]
assert sfnomdata["key"]=="sf"
assert sfupdata["key"]=="sfup"
assert sfdowndata["key"]=="sfdown"

# initialize structure
npdata = {'ptbelow20': {}, 'ptabove20':{}}

# loop over nominal, up and down
for systag, sfdata in zip(['nominal', 'up', 'down'], [sfnomdata, sfupdata, sfdowndata]):

    # further select relevant part of the data
    sfdata = sfdata["value"]["content"]
    sfbelow20 = sfdata[0]
    sfabove20 = sfdata[1]
    assert sfbelow20["key"]=="RecoBelow20"
    assert sfabove20["key"]=="RecoAbove20"
    
    # loop over pt bins
    for pttag, sfdata in zip(['ptbelow20', 'ptabove20'], [sfbelow20, sfabove20]):
        sfdata = sfdata["value"]
        eta_edges = np.array(sfdata["edges"][0])
        pt_edges = np.array(sfdata["edges"][1])
        values = np.array(sfdata["content"])
        values = np.reshape(values, (len(eta_edges)-1,len(pt_edges)-1))
        npdata[pttag][systag] = (eta_edges, pt_edges, values)
      
# loop over pt bins
for pttag in ['ptbelow20', 'ptabove20']:
    # convert to TH2
    eta_edges = npdata[pttag]['nominal'][0]
    eta_edges = np.clip(eta_edges, -2.5, 2.5)
    pt_edges = npdata[pttag]['nominal'][1]
    pt_edges = np.clip(pt_edges, 0., 200.)
    nomvalues = npdata[pttag]['nominal'][2]
    upvalues = npdata[pttag]['up'][2]
    downvalues = npdata[pttag]['down'][2]
    sys = np.maximum(abs(upvalues-nomvalues), 
                     abs(downvalues-nomvalues))
    hist = ROOT.TH2D(pttag, 'Electron reco scale factors',
                     len(eta_edges)-1, array.array('f', list(eta_edges)),
                     len(pt_edges)-1, array.array('f', list(pt_edges)))
    for i in range(1,len(eta_edges)):
        for j in range(1,len(pt_edges)):
            hist.SetBinContent(i,j, nomvalues[i-1,j-1])
            hist.SetBinError(i,j, sys[i-1,j-1])

    # make a plot
    figname = args.sffile.replace('.json','_reco_{}.png'.format(pttag))
    ptlabels = {'ptbelow20': 'pT < 20 GeV', 'ptabove20': 'pT > 20 GeV'}
    title = 'Electron reco scale factor for {}, {}'.format(args.year, ptlabels[pttag])
    plot2dhistogram( hist, figname, histtitle=title, cmin=0.9, cmax=1.1 ) 
