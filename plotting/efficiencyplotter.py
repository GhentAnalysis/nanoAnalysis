#############################################################################
# Plot data and MC efficiencies (e.g. trigger efficiencies) and their ratio #
#############################################################################

import sys
import os
import ROOT
import numpy as np
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
import plotting.plottools as pt
import plotting.histplotter as hp

def drawextrainfo(infos):
    # write extra info to current plot
    # infos is a list of dicts
    # minimal dict entries are 'text','posx' and 'posy'
    # additional optional entries: see below
    tinfo = ROOT.TLatex()
    for info in infos:
      keys = info.keys()
      if 'textcolor' in keys: 
        tinfo.SetTextColor(info['textcolor'])
      if 'textsize' in keys: 
        tinfo.SetTextSize(info['textsize'])
      if 'textfont' in keys: 
        tinfo.SetTextFont(info['textfont'])
      tinfo.DrawLatexNDC(info['posx'],info['posy'],info['text'])

def getweightedvalue(valuegraph, spectrumhist):
    # average out the per-bin values and uncertainties in valuegraph, 
    # using the spectrumhist as weighting factors.
    # make sure that the bins in spectrumhist correspond to the points in valuegraph!
    
    avgval = avgup = avgdown = 0
    sumweighting = 0.
    for i in range(valuegraph.GetN()):
        # check that bin and point align
        if valuegraph.GetX()[i] != spectrumhist.GetBinCenter(i+1):
            raise Exception('ERROR in efficienceyplotter.getweightedvalue:'
                  +' points in graph and bin centers in hist do not correspond')
        # determine the weight factor
        weightfactor = spectrumhist.GetBinContent(i+1)
        sumweighting += weightfactor
        # determine the value for this bin
        avgval += valuegraph.GetY()[i]*weightfactor
        avgup += np.power(valuegraph.GetErrorYhigh(i),2)*weightfactor
        avgdown += np.power(valuegraph.GetErrorYlow(i),2)*weightfactor
        # printouts for testing
        #print('--- point '+str(i)+' ---')
        #print('weightfactor: '+str(weightfactor))
        #print('value: '+str(valuegraph.GetY()[i]))
        #print('uperror:'+str(valuegraph.GetErrorYhigh(i)))
        #print('uperror:'+str(valuegraph.GetErrorYlow(i)))
    # calculate the averages
    avgval /= sumweighting
    avgup = np.sqrt(avgup/sumweighting)
    avgdown = np.sqrt(avgdown/sumweighting)
    return (avgval,avgup,avgdown)

def plotefficiencies(datagraphlist, figname, 
                    simgraph=None, datahistlist=None, simhist=None, simsysthist=None,
                    datacolorlist=None, datalabellist=None, simcolor=None, simlabel=None,
                    yaxtitle=None, yaxtitlesize=None, yaxtitleoffset=None,
                    xaxtitle=None, xaxtitlesize=None, xaxtitleoffset=None,
                    extracmstext='', lumi=None,
                    dosimavg=False, dodataavg=False):
    ### plot efficiencies
    # input arguments:
    # - datagraphlist: list of TGraphAsymmError containing data efficiencies
    #   note: the last element in datagraphlist will be highlighted,
    #         i.e. thicker line and bigger marker;
    #         also average value (if dodataavg is True) is calculated for the last one.
    # - simgraph: TGraphAsymmError containing simulated efficiencies.
    # - datahistlist: list of TH1D containing absolute data yields instead of efficiencies.
    #        if datahistlist is not None it is assumed to match datagraphlist,
    #   but individual elements can still be None in order not to plot them.
    # - simhist: TH1D containing absolute simulation yield instead of efficiency.
    # - simsyshist: TH1D containing systematic variations to be applied to simgraph.
    # others: plot aesthetics.
    
    pt.setTDRstyle()
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    ### define global parameters for size and positioning
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    cheight = 600 # height of canvas
    cwidth = 600 # width of canvas
    rfrac = 0.25 # fraction of ratio plot in canvas
    # fonts and sizes:
    labelfont = 4; labelsize = 22
    axtitlefont = 4
    if xaxtitlesize is None: xaxtitlesize = 22
    if yaxtitlesize is None: yaxtitlesize = 22
    infofont = 4; infosize = 20
    # title offset
    if yaxtitleoffset is None: yaxtitleoffset = 1.9
    if xaxtitleoffset is None: xaxtitleoffset = 4.5
    # margins:
    p1topmargin = 0.07
    p2bottommargin = 0.4
    leftmargin = 0.15
    rightmargin = 0.05
    # legend box
    p1legendbox = [leftmargin+0.03,1-p1topmargin-0.25,1-rightmargin-0.03,1-p1topmargin-0.03]
    p2legendbox = [leftmargin+0.03,0.84,1-rightmargin-0.03,0.97]
    # marker properties for data
    markerstyle = 20
    markersize = 0.5
    # y axis range
    yaxmin = 0. # note: spectrum plotting will not work well for yaxmin > 0
                # could introduce offset but this gives issues with weighting in averaging
    yaxmax = 1.5
    # x axis range
    xaxmin = simgraph.GetX()[0]-simgraph.GetErrorXlow(0)
    xaxmax = simgraph.GetX()[simgraph.GetN()-1]+simgraph.GetErrorXhigh(simgraph.GetN()-1)
    # default colors
    if datacolorlist is None: datacolorlist = [ROOT.kBlue]*len(datahistlist)
    if simcolor is None: simcolor = ROOT.kBlack

    ### set properties of absolute simulation histogram showing the spectrum
    csurface = (yaxmax-yaxmin)*(xaxmax-xaxmin)
    if simhist is not None:
        simhist.SetFillColor(ROOT.kGray)
        simhist.SetLineColor(simcolor)
        simhist.SetLineWidth(1)
        # normalization and offset
        integral = simhist.Integral("width")
        scale = 0.15*csurface/integral
        simhist.Scale(scale)
        for i in range(1,simhist.GetNbinsX()+1):
            # set bins with zero error to zero.
            # (they could be filled with arbitrary values
            #  to avoid effiency errors but dont plot these)
            if simhist.GetBinError(i)==0:
                simhist.SetBinContent(i,0)

    ### set properties of absolute data histograms showing the spectrum
    if datahistlist is not None:
        for i,datahist in enumerate(datahistlist):
            if datahist is not None:
                color = datacolorlist[i]
                datahist.SetLineWidth(1)
                datahist.SetLineColor(color)
                # highlight the last graph
                if i==len(datahistlist)-1:
                    datahist.SetLineWidth(2)
                # normalization and offset
                integral = datahist.Integral("width")
                scale = 0.15*csurface/integral
                datahist.Scale(scale)
                for i in range(1,datahist.GetNbinsX()+1):
                    # set bins with zero error to zero.
                    # (they could be filled with arbitrary values
                    #  to avoid effiency errors but dont plot these)
                    if datahist.GetBinError(i)==0:
                        datahist.SetBinContent(i,0) 

    ### calculate statistical and total simulation error and set graph properties
    # note: simstaterror is a copy of simgraph before modifying it, containing statistical errors;
    #       simerror starts as a copy of simgraph but contains both statistical and systematic errors;
    #            simgraph gets is vertical errors set to zero (for plotting purposes)
    simstaterror = simgraph.Clone()
    simerror = simgraph.Clone()
    # if a systematics histogram is given, calculate total stat+sys error and modify simerror
    if simsysthist is not None:
        for i in range(0,simgraph.GetN()):
            staterrorup = simgraph.GetErrorYhigh(i)
            staterrordown = simgraph.GetErrorYlow(i)
            systerror = mcsyshist.GetBinContent(i+1)
            toterrorup = np.sqrt(np.power(staterrorup,2)+np.power(systerror,2))
            toterrordown = np.sqrt(np.power(staterrordown,2)+np.power(systerror,2))
            toterrorup = min(1-simgraph.GetY()[i],toterrorup)
            toterrordown = min(simgraph.GetY()[i],toterrordown)
            simerror.SetPointEYhigh(i,toterrorup)
            simerror.SetPointEYlow(i,toterrordown)
    # in any case, put vertical errors of simgraph to zero
    for i in range(0,simgraph.GetN()):
        simgraph.SetPointEYhigh(i,0.)
        simgraph.SetPointEYlow(i,0.)
    # set drawing properties of simgraph and simerror (simstaterror is not drawn, maybe introduce later)
    simgraph.SetLineColor(ROOT.kBlack)
    simgraph.SetLineWidth(2)
    simerror.SetFillStyle(3244)
    simerror.SetLineWidth(0)
    simerror.SetFillColor(ROOT.kGray+2)
    simerror.SetMarkerStyle(0)

    ### calculate statistical and total mc error (scaled to unit bin content) 
    scstaterror = simstaterror.Clone()
    scerror = simerror.Clone()
    for i in range(0,simgraph.GetN()):
        scstaterror.SetPoint(i,simgraph.GetX()[i],1.)
        scerror.SetPoint(i,simgraph.GetX()[i],1.)
        if not simgraph.GetY()[i]==0:
            scstaterror.SetPointEYhigh(i,simstaterror.GetErrorYhigh(i)/simstaterror.GetY()[i])
            scstaterror.SetPointEYlow(i,simstaterror.GetErrorYlow(i)/simstaterror.GetY()[i])
            scerror.SetPointEYhigh(i,simerror.GetErrorYhigh(i)/simerror.GetY()[i])
            scerror.SetPointEYlow(i,simerror.GetErrorYlow(i)/simerror.GetY()[i])
        else:
            scstaterror.SetPointEYhigh(i,0.)
            scstaterror.SetPointEYlow(i,0.)
            scerror.SetPointEYhigh(i,0.)
            scerror.SetPointEYlow(i,0.)
    scstaterror.SetFillStyle(1001)
    scerror.SetFillStyle(3254)
    scstaterror.SetFillColor(ROOT.kGray+1)
    scerror.SetFillColor(ROOT.kBlack)
    scstaterror.SetMarkerStyle(1)
    scerror.SetMarkerStyle(1)

    ### operations on data histogram
    for i,datagraph in enumerate(datagraphlist):
        color = datacolorlist[i]
        datagraph.SetMarkerStyle(markerstyle)
        datagraph.SetMarkerColor(color)
        datagraph.SetMarkerSize(markersize)
        datagraph.SetLineWidth(1)
        datagraph.SetLineColor(color)
        # highlight the last graph
        if i==len(datagraphlist)-1:
            datagraph.SetLineWidth(2)
            datagraph.SetMarkerSize(markersize*2)

    ### calculate data to mc ratio
    ratiographlist = []
    for datagraph in datagraphlist:
        ratiograph = datagraph.Clone()
        for i in range(0,simgraph.GetN()):
            if not simgraph.GetY()[i]==0:
                ratiograph.GetY()[i] *= 1./simgraph.GetY()[i]
                ratiograph.SetPointEYhigh(i,datagraph.GetErrorYhigh(i)/simgraph.GetY()[i])
                ratiograph.SetPointEYlow(i,datagraph.GetErrorYlow(i)/simgraph.GetY()[i])
            # avoid drawing empty mc or data bins
            else: ratiograph.GetY()[i] = 1e6
            if(datagraph.GetY()[i]<=0): ratiograph.GetY()[i] += 1e6
        ratiographlist.append(ratiograph)

    ### get weighted efficiency values
    info = []
    if( dosimavg and simhist is not None ):
        mctup = getweightedvalue(simerror,simhist)
        info.append( { 'text':'weighted sim eff:',
                        'posx':0.65,'posy':0.5,
                        'textsize':infosize*0.75,'textfont':infofont} )
        info.append( { 'text':'{:.4}'.format(mctup[0])
                        +'^{+'+'{:.4}'.format(mctup[1])+'}'
                        +'_{-'+'{:.4}'.format(mctup[2])+'}',
                        'posx':0.65,'posy':0.45,
                        'textsize':infosize,'textfont':infofont} )
    if( dodataavg and datahistlist is not None ):
        datahist = datahistlist[-1]
        datagraph = datagraphlist[-1]
        datatup = getweightedvalue(datagraph,datahist)
        info.append( { 'text':'weighted data eff:',
                        'posx':0.65,'posy':0.35,
                        'textsize':infosize*0.75,'textfont':infofont,
                        'textcolor':ROOT.kRed} )
        info.append( { 'text':'{:.4}'.format(datatup[0])
                        +'^{+'+'{:.4}'.format(datatup[1])+'}'
                        +'_{-'+'{:.4}'.format(datatup[2])+'}',
                        'posx':0.65,'posy':0.3,
                        'textsize':infosize,'textfont':infofont,
                        'textcolor':ROOT.kRed} )
                
    ### make legend for upper plot and add all histograms
    legend = ROOT.TLegend(p1legendbox[0],p1legendbox[1],p1legendbox[2],p1legendbox[3])
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    if simlabel is not None: legend.AddEntry(simgraph,simlabel,"l")
    legend.AddEntry(simerror,"Total sim. unc.","f")
    if datalabellist is not None:
        for i,datagraph in enumerate(datagraphlist):
            legend.AddEntry(datagraph,datalabellist[i],"lpe")
    if( simlabel is not None and simhist is not None ):
        legend.AddEntry(simhist,simlabel)

    ### make legend for lower plot and add all histograms
    legend2 = ROOT.TLegend(p2legendbox[0],p2legendbox[1],p2legendbox[2],p2legendbox[3])
    legend2.SetNColumns(2); 
    legend2.SetFillStyle(0);
    legend2.AddEntry(scstaterror, "Stat. uncertainty", "f");
    legend2.AddEntry(scerror, "Uncertainty", "f");

    ### make canvas and pads
    c1 = ROOT.TCanvas("c1","c1")
    c1.SetCanvasSize(cwidth,cheight)
    pad1 = ROOT.TPad("pad1","",0.,rfrac,1.,1.)
    pad1.SetTopMargin(p1topmargin)
    pad1.SetBottomMargin(0.03)
    pad1.SetLeftMargin(leftmargin)
    pad1.SetRightMargin(rightmargin)
    pad1.SetGrid()
    pad1.Draw()
    pad2 = ROOT.TPad("pad2","",0.,0.,1.,rfrac)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(p2bottommargin)
    pad2.SetLeftMargin(leftmargin)
    pad2.SetRightMargin(rightmargin)
    pad2.SetGrid()
    pad2.Draw()
    
    ### make upper part of the plot
    pad1.cd()
    simerror.SetMinimum(yaxmin)
    simerror.SetMaximum(yaxmax)

    # X-axis layout
    xax = simerror.GetXaxis()
    xax.SetNdivisions(5,4,0,ROOT.kTRUE)
    xax.SetLimits(xaxmin,xaxmax)
    xax.SetLabelSize(0)
    # Y-axis layout
    yax = simerror.GetYaxis()
    yax.SetMaxDigits(3)
    yax.SetNdivisions(8,4,0,ROOT.kTRUE)
    yax.SetLabelFont(10*labelfont+3)
    yax.SetLabelSize(labelsize)
    if yaxtitle is not None:
        yax.SetTitle(yaxtitle)
        yax.SetTitleFont(10*axtitlefont+3)
        yax.SetTitleSize(yaxtitlesize)
        yax.SetTitleOffset(yaxtitleoffset)

    # draw simerror first to get range correct
    simerror.Draw("A 2") # option A to draw axes, option 2 to draw error as rectangular bands
    if simhist is not None: simhist.Draw("hist same")
    if datahistlist is not None:
        for datahist in datahistlist:
            if datahist is not None: datahist.Draw("hist same")
    # draw data graphs
    # note: drawing two times to have perpendicular bars at vertical but not at horizontal error bars
    for datagraph in datagraphlist: 
        datagraph.Draw("Z P same") # option Z to suppress perpendicular bars
        for i in range(datagraph.GetN()): 
            datagraph.SetPointEXhigh(i,0)
            datagraph.SetPointEXlow(i,0)
        datagraph.Draw("P same")
    simgraph.Draw("Z same") # option Z to suppress perpendicular bars
    legend.Draw("same")
    ROOT.gPad.RedrawAxis()

    # make and draw unit line
    line = ROOT.TLine(xaxmin,1,xaxmax,1)
    line.SetLineStyle(2)
    line.Draw("same")

    # draw extra info
    if len(info)>0:
        drawextrainfo(info)

    # draw header
    lumistr = ''
    if lumi is not None:
        lumistr = '{0:.3g}'.format(lumi/1000.)+' fb^{-1} (13 TeV)'
    pt.drawLumi(pad1,cms_in_grid=False,
                extratext=extracmstext,lumitext=lumistr)

    ### make the lower part of the plot
    pad2.cd()
    # X-axis layout
    xax = scerror.GetXaxis()
    xax.SetLimits(xaxmin,xaxmax)
    xax.SetNdivisions(5,4,0,ROOT.kTRUE)
    xax.SetLabelSize(labelsize)
    xax.SetLabelFont(10*labelfont+3)
    if xaxtitle is not None:
        xax.SetTitle(xaxtitle)
        xax.SetTitleFont(10*axtitlefont+3)
        xax.SetTitleSize(xaxtitlesize)
        xax.SetTitleOffset(xaxtitleoffset)
    # Y-axis layout
    yax = scerror.GetYaxis()
    #yax.SetLimits(0.9,1.0999)
    scerror.SetMinimum(0.95)
    scerror.SetMaximum(1.0499)
    yax.SetMaxDigits(3)
    yax.SetNdivisions(4,5,0)
    yax.SetLabelFont(10*labelfont+3)
    yax.SetLabelSize(labelsize)
    yax.SetTitle("Obs./pred.")
    yax.SetTitleFont(10*axtitlefont+3)
    yax.SetTitleSize(yaxtitlesize)
    yax.SetTitleOffset(yaxtitleoffset)

    # draw objects
    scerror.Draw("A 2")
    scstaterror.Draw("2 same")
    # draw ratio graphs
    # note: drawing two times to have perpendicular bars at vertical but not at horizontal error bars
    for ratiograph in ratiographlist:
        ratiograph.Draw("Z P same") # option Z to suppress perpendicular bars
        for i in range(ratiograph.GetN()):
            ratiograph.SetPointEXhigh(i,0)
            ratiograph.SetPointEXlow(i,0)
        ratiograph.Draw("P same")
    legend2.Draw("same")
    ROOT.gPad.RedrawAxis()

    # make and draw unit ratio line
    line2 = ROOT.TLine(xaxmin,1,xaxmax,1)
    line2.SetLineStyle(2)
    line2.Draw("same")
    
    ### save the plot
    c1.SaveAs(figname)
