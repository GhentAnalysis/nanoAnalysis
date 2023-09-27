########################################################################
# some small tools for working with histograms and lists of histograms #
########################################################################

# attempt to fully switch to uproot instead of PyROOT

# import python modules
import sys
import os
import numpy as np
import uproot
# import local modules
from pathlib import Path
sys.path.append(str(Path(__file__).parents[2]))
import tools.listtools as lt
from tools.histogram import Histogram


### histogram name reading ###

def loadhistnames(histfile,
                  mustcontainall=[], mustcontainone=[],
                  maynotcontainall=[],maynotcontainone=[],
                  allow_tgraphs=False):
    ### read a root file containing histograms and make a list of histogram names.
    # note: objects are not loaded (for speed), only a list of names is retrieved.
    with uproot.open(histfile) as f:
        classnames = f.classnames()
    histnames = ([key for key,val in classnames.items()
                  if val.startswith('TH')])
    if allow_tgraphs:
        graphnames = ([key for key,val in classnames.items()
                       if val.startswith('TGraph')])
        histnames += graphnames
    histnames = lt.subselect_strings(histnames,
      mustcontainall=mustcontainall, mustcontainone=mustcontainone,
      maynotcontainall=maynotcontainall, maynotcontainone=maynotcontainone)[1]
    return histnames

def loadallhistnames(histfile, allow_tgraphs=False):
    ### alias for loadhistnames with no selection applied
    return loadhistnames(histfile, allow_tgraphs=allow_tgraphs)


### histogram reading ###

def loadhistogramlist(histfile, histnames, do_checks=True, reset_names=False):
    ### load histograms specified by name from a file
    hists = []
    with uproot.open(histfile) as f:
        classnames = f.classnames()
        for histname in histnames:
            if do_checks:
                if histname not in classnames.keys():
                    msg = 'WARNING in loadhistogramlist:'
                    msg += ' name {} not found in input file {}.'.format(histname, histfile)
                    print(msg)
                    continue
                if not classnames[histname].startswith('TH'):
                    msg = 'WARNING in loadhistogramlist:' 
                    msg += ' name {} does not seem to be a THx.'.format(histname)
                    print(msg)
            hist = Histogram.from_uproot(f[histname])
            if reset_names: hist.name = histname
            hists.append(hist)
    return hists

def loadallhistograms(histfile, reset_names=False, allow_tgraphs=False):
    ### read a root file containing histograms and load all histograms to a list
    histnames = loadallhistnames(histfile, allow_tgraphs=allow_tgraphs)
    return loadhistogramlist(histfile, histnames, do_checks=False, reset_names=reset_names)

def loadhistograms(histfile, 
                   mustcontainall=[], mustcontainone=[],
                   maynotcontainall=[], maynotcontainone=[],
                   reset_names=False, allow_tgraphs=False):
    ### read a root file containing histograms and load selected histograms to a list.
    histnames = loadhistnames(histfile,
      mustcontainone=mustcontainone, mustcontainall=mustcontainall,
      maynotcontainone=maynotcontainone, maynotcontainall=maynotcontainall,
      allow_tgraphs=allow_tgraphs)
    return loadhistogramlist(histfile, histnames, do_checks=False, reset_names=reset_names)


### histogram subselection ###

def selecthistograms(histlist,
                     mustcontainone=[], mustcontainall=[],
                     maynotcontainone=[], maynotcontainall=[]):
    idlist = [hist.name for hist in histlist]
    (indlist,selhistlist) = lt.subselect_objects(histlist,idlist,
        mustcontainone=mustcontainone,mustcontainall=mustcontainall,
        maynotcontainone=maynotcontainone,maynotcontainall=maynotcontainall)
    return (indlist,selhistlist)

def findhistogram(histlist, name):
    ### find a histogram with a given name in a list
    # returns None if no histogram with the requested name is found
    for hist in histlist:
        if hist.name==name: return hist
    return None


### histogram clipping ###

def cliphistogram(hist, clipboundary=0):
    ### clip a histogram to minimum zero
    # also allow a clipboundary different from zero, useful for plotting 
    # (e.g. to ignore artificial small values)
    values = hist.values(flow=True)
    errors = hist.errors(flow=True)
    inds = np.nonzero(values < clipboundary)
    values[inds] = 0.
    errors[inds] = 0.
    if np.sum(values) < 1e-12: values[1] = 1e-6
    hist.values = values
    hist.errors = errors

def cliphistograms(histlist, clipboundary=0):
    ### apply cliphistogram on all histograms in a list
    for hist in histlist: cliphistogram(hist, clipboundary=clipboundary)

def clipallhistograms(histfile, mustcontainall=[], clipboundary=0):
    ### apply cliphistogram on all histograms in a file
    histlist = loadallhistograms(histfile)
    if len(mustcontainall)==0:
        cliphistograms(histlist, clipboundary=clipboundary)
    else:
        (indlist,_) = selecthistograms(histlist, mustcontainall=mustcontainall)
        for index in indlist: cliphistogram(histlist[index], clipboundary=clipboundary)
    tempfilename = histfile[:-5]+'_temp.root'
    with uproot.recreate(tempfilename) as f:
        # writing of custom Histogram objects not yet implemented,
        # the line below will likely fail or give nonsense
        for hist in histlist: f[hist.name] = hist
    os.system('mv '+tempfilename+' '+histfile)


### finding minimum and maximum ###

def getminmax(histlist, includebinerror=False):
    # get suitable minimum and maximum values for plotting a hist collection (not stacked)
    allvalues = []
    for hist in histlist:
        if includebinerror:
          allvalues.append(hist.values() + hist.errors())
          allvalues.append(hist.values() - hist.errors())
        else: allvalues.append(hist.values())
    allvalues = np.concatenate(tuple(allvalues))
    return (np.min(allvalues), np.max(allvalues))

def getminmaxmargin(histlist, includebinerror=False, clip=False):
    (totmin,totmax) = getminmax(histlist, includebinerror=includebinerror)
    topmargin = (totmax-totmin)/2.
    bottommargin = (totmax-totmin)/5.
    minv = totmin-bottommargin
    maxv = totmax+topmargin
    if( clip and minv<0 ): minv = 0
    return (minv,maxv)


### histogram conversion ###

def tgraphtohist( graph ):
    raise Exception('Not yet implemented')

    # get list of x values and sort them
    xvals = []
    for i in range(graph.GetN()): xvals.append(graph.GetX()[i])
    xvals = np.array(xvals)
    sortedindices = np.argsort(xvals)
    # make bins
    bins = []
    for i in sortedindices: bins.append(graph.GetX()[i]-graph.GetErrorXlow(i))
    bins.append(graph.GetX()[i]+graph.GetErrorXhigh(i))
    # make histogram
    hist = ROOT.TH1D("","",len(bins)-1,array('f',bins))
    # set bin contents
    for i in range(1,hist.GetNbinsX()+1):
        bincontent = graph.GetY()[sortedindices[i-1]]
        binerror = max(graph.GetErrorYlow(sortedindices[i-1]),
                        graph.GetErrorYhigh(sortedindices[i-1]))
        hist.SetBinContent(i,bincontent)
        hist.SetBinError(i,binerror)
    hist.name = graph.name
    hist.title = graph.title
    return hist


### histogram calculations ###

def binperbinmaxvar( histlist, nominalhist ):
    ### get the bin-per-bin maximum absolute variation of histograms w.r.t. a nominal histogram
    varvalues = np.zeros((len(histlist),len(nominalhist.values)))
    for i,hist in enumerate(histlist):
        varvalues[i,:] = abs(hist.values() - nominalhist.values())
    maxvarvalues = np.max(varvalues, axis=0)
    return maxvarvalues # or return Histogram with these values?

def envelope( histlist, returntype='tuple' ):
    ### return two histograms that form the envelope of all histograms in histlist.
    # arguments:
    # - returntype: if 'tuple', returns tuple of two histograms (lower bound and upper bound);
    #               if 'hist', returns a single histogram with same lower and upper bounds
    #               (and bin contents chosen symmetrically between them).
    if( len(histlist)<2 ):
        msg = 'ERROR in histtools.envelope: at least two histograms required.'
        raise Exception(msg)
    allvalues = np.zeros((len(histlist),len(histlist[0].values)))
    for i,hist in enumerate(histlist):
        allvalues[i,:] = hist.values()
    lower = np.min(allvalues, axis=0)
    upper = np.max(allvalues, axis=0)
    if returntype=='tuple': return (lower, upper)
    elif returntype=='hist':
        bounds = np.vstack((lower, upper))
        values = np.mean(bounds, axis=0)
        errors = upper - values
        return Histogram(values=values, errors=errors)
    else:
        msg = 'ERROR in histtools.envelope:'
        msg += ' return type {} not recognized.'.format(returntype)
        raise Exception(msg)
    
def rootsumsquare( histlist ):
    ### return a histogram that is the root-sum-square of all histograms in histlist.
    # check the input list
    if( len(histlist)<1 ):
        msg = 'ERROR in histtools.rootsumsquare: at least one histogram required.'
        raise Exception(msg)
    allvalues = np.zeros((len(histlist),len(histlist[0].values)))
    for i,hist in enumerate(histlist):
        allvalues[i,:] = hist.values()
    rss = np.sqrt(np.sum(np.square(allvalues), axis=0))
    return rss
