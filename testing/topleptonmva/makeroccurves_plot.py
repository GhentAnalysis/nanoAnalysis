################################################
# Script for making ROC curves for lepton MVAs #
################################################
# part 2: reading histograms from file and plotting

# imports
import sys
import os
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt

# local imports
sys.path.append(str(Path(__file__).parents[2]))
import tools.histtools as ht


def getefffromhist( hist ):
    ### internal helper function for calculating ROC efficienies
    # input arguments:
    # - hist: TH1 containing mva scores
    # returns:
    # numpy array with efficiencies (in increasing order)
    nbins = hist.GetNbinsX()
    xlow = hist.GetBinLowEdge(1)
    xhigh = hist.GetBinLowEdge(nbins) + hist.GetBinWidth(nbins)
    eff = np.zeros(nbins+1)
    # loop over bins
    for i in range(0,nbins+1):
        # take integral starting from next bin to end
        # (note that the range in TH1.Integral is inclusive)
        eff[i] = hist.Integral(i+1, nbins+1)
    # normalize 
    eff = eff[::-1]/eff[0]
    # return result
    return eff

def getrocfromhists( signalhist, backgroundhist ):
    ### internal helper function for calculating ROC curves
    # input arguments:
    # - signalhist: TH1 containing mva scores for signal
    # - backgroundhist: TH1 containing mva scores for background
    # returns:
    # a dict with following fields:
    # - effs: numpy array with signal efficiencies (in increasing order)
    # - effb: numpy array with background efficiencies (in increasing order)
    effs = getefffromhist(signalhist)
    effb = getefffromhist(backgroundhist)
    if( len(effs)!=len(effb) ):
        msg = 'ERROR in getrocfromhist:'
        msg += ' number of bins in signal and background do not agree.'
        raise Exception(msg)
    return {'effs': effs, 'effb': effb}

def plotrocs( inputlist, title=None,
    xaxtitle='Background efficiency',
    yaxtitle='Signal efficiency',
    xaxlog=False, yaxlog=False,
    xaxrange=None, yaxrange=None ):
    ### make ROC plots
    # input arguments:
    # - inputlist: list containing dicts (one dict per curve)
    #   with entries 'effs','effb','color','label'
    fig,ax = plt.subplots()
    for inputdict in inputlist:
        ax.plot(inputdict['effb'],inputdict['effs'],
        label=inputdict['label'],color=inputdict['color'],
        linewidth=2.)
    ax.legend(loc='lower right')
    if title is not None: ax.set_title(title, fontsize=15)
    # general axis titles:
    ax.set_xlabel(xaxtitle, fontsize=15)
    ax.set_ylabel(yaxtitle, fontsize=15)
    if xaxlog: ax.set_xscale('log')
    if yaxlog: ax.set_yscale('log')
    if xaxrange is not None: ax.set_xlim(*xaxrange)
    if yaxrange is not None: ax.set_ylim(*yaxrange)
    ax.grid()
    return (fig,ax)


if __name__=='__main__':

    # input arguments:
    parser = argparse.ArgumentParser(description='Make ROC curves for lepton MVAs')
    parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
    parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # make output directory
    if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

    # load all histograms from input file
    # note: histograms are assumed to be named flavour_ptLowtoHigh_mvaVar_prompt/nonprompt
    hists = ht.loadallhistograms(args.inputfile)
    histdict = {}
    histnames = []
    for hist in hists:
        name = hist.GetName()
        histdict[name] = hist
        histnames.append(name)
    phistnames = sorted([el for el in histnames if el.endswith('_prompt')])
    nphistnames = sorted([el for el in histnames if el.endswith('_nonprompt')])

    # find which instances are present
    elements = [el.split('_') for el in phistnames]
    flavours = list(set([el[0] for el in elements]))
    ptcuts = list(set([el[1] for el in elements]))
    mvavars = list(set([el[2] for el in elements]))
    print('Found following instances in input file:')
    print('- flavours: {}'.format(flavours))
    print('- pt cuts: {}'.format(ptcuts))
    print('- MVAs: {}'.format(mvavars))

    # define colors and labels
    colors = {
      'mvaTTH': 'b',
      'mvaTOP': 'r'
    }
    labels = {
      'mvaTTH': 'mvaTTH',
      'mvaTOP': 'mvaTOP'
    }
    flavourlabels = {
        'muon': 'Muons',
        'electron': 'Electrons'
    }
    ptlabels = {
        'pt10toInf': 'pT > 10 GeV',
        'pt10to25': '10 < pT < 25 GeV',
        'pt25toInf': 'pT > 25 GeV'
    }

    # loop over instances
    for flavour in flavours:
        for ptcut in ptcuts:
            mvarocs = []
            for mvavar in mvavars:
                print('Making ROC curve for {}, {}, {}'.format(flavour,ptcut,mvavar))
                # get the correct histograms
                name = '_'.join([flavour, ptcut, mvavar])
                phist = histdict.get(name+'_prompt', None)
                nphist = histdict.get(name+'_nonprompt', None)
                if( phist is None or nphist is None ):
                    msg = 'WARNING: expeced histograms not found,'
                    msg += ' skipping this instance.'
                    print(msg)
                    continue
                roc = getrocfromhists(phist, nphist)
                mvarocs.append({
                    'effs':roc['effs'],
                    'effb':roc['effb'],
                    'label':labels[mvavar],
                    'color':colors[mvavar]})
            # make the figure
            fig,ax = plotrocs( mvarocs, xaxlog=True,
                       xaxrange=(1e-3,1), yaxrange=(0.7,1.001) )
            # additional figure editing
            ax.text(0.96, 0.25, flavourlabels.get(flavour,flavour),
                    ha='right', transform=ax.transAxes)
            ax.text(0.96, 0.2, ptlabels.get(ptcut,ptcut),
                    ha='right', transform=ax.transAxes)
            figname = '_'.join([flavour,ptcut])+'.png'
            figname = os.path.join(args.outputdir, figname)
            fig.savefig(figname)
