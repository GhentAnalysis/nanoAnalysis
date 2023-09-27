#########################################################################
# plot and compare multiple theory distributions to a data distribution #
#########################################################################

import ROOT
import sys
import os
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
import plotting.plottools as pt
import tools.histtools as ht

def plotdifferential(
    theoryhists,
    datahist,
    systhists=None,
    statdatahist=None,
    figname=None, title=None, xaxtitle=None, yaxtitle=None,
    dolegend=True, labellist=None, 
    colorlist=None,
    logy=False, ymaxlinfactor=1.8, yminlogfactor=0.2, ymaxlogfactor=100,
    drawoptions='', 
    lumitext='', extracmstext = '',
    ratiorange=None, ylims=None, yminzero=False,
    extrainfos=[], infosize=None, infoleft=None, infotop=None ):
    ### plot multiple overlaying histograms (e.g. for shape comparison)
    # note: the ratio plot will show ratios w.r.t. the first histogram in the list!
    # arguments:
    # - theoryhists, colorlist, labellist: lists of TH1, ROOT colors and labels respectively
    # - systhists: list of TH1 with systematic uncertainties on theoryhists
    # - datahist: TH1 with the measurement values
    # - statdatahist: TH1 with same bin contents as datahist but statistical-only errors
    # - figname: name of the figure to save (if None, do not save but return plot dictionary)
    # - title, xaxtitle, yaxtitle, figname: self-explanatory
    # - dolegend: boolean whether to make a legend (histogram title is used if no labellist)
    # - logy: boolean whether to make y-axis logarithmic
    # - ymaxlinfactor: factor by which to multiply maximum y value (for linear y-axis)
    # - yminlogfactor and ymaxlogfactor: same as above but for log scale
    # - drawoptions: string passed to TH1.Draw
    #   see https://root.cern/doc/master/classTHistPainter.html for a full list of options
    # - lumitext and extracmstext: luminosity value and extra text
    # - ratiorange: a tuple of (ylow,yhigh) for the ratio pad, default (0,2)
    # - ylims: a tuple of (ylow,yhigh) for the upper pad
    # - yminzero: whether to clip minimum y to zero.
    # - extrainfos is a list of strings with extra info to display
    # - infosize: font size of extra info
    # - infoleft: left border of extra info text (default leftmargin + 0.05)
    # - infotop: top border of extra info text (default 1 - topmargin - 0.1)

    pt.setTDRstyle()
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    ### parse arguments
    if colorlist is None:
        colorlist = ([ROOT.kAzure-4, ROOT.kAzure+6, ROOT.kViolet, ROOT.kMagenta-9,
                      ROOT.kRed, ROOT.kGreen+1, ROOT.kGreen-1])
    if( len(theoryhists)>len(colorlist) ):
        raise Exception('ERROR in plotdifferential:'
            +' theory histogram list is longer than color list')
    if(labellist is not None and len(labellist)!=len(theoryhists)):
        raise Exception('ERROR in plotdifferential:'
            +' length of label list does not agree with theory histogram list')
    if( systhists is not None and len(systhists)!=len(theoryhists) ):
        raise Exception('ERROR in plotdifferential:'
            +' length of systematic histogram list does not agree with theory histogram list')

    ### define global parameters for size and positioning
    cheight = 600 # height of canvas
    cwidth = 600 # width of canvas
    rfrac = 0.25 # fraction of bottom plot showing the ratio
    # fonts and sizes:
    labelfont = 4; labelsize = 22
    axtitlefont = 4; axtitlesize = 26
    infofont = 4
    if infosize is None: infosize = 20
    legendfont = 4
    # margins and title offsets
    ytitleoffset = 1.5
    p1topmargin = 0.07 
    p1bottommargin = 0.03
    xtitleoffset = 3.5
    p2topmargin = 0.01
    p2bottommargin = 0.4
    leftmargin = 0.15
    rightmargin = 0.05
    # legend box
    pentryheight = 0.07
    nentries = 1 + len(theoryhists)
    if nentries>3: pentryheight = pentryheight*0.8
    plegendbox = ([leftmargin+0.5,1-p1topmargin-0.03-pentryheight*nentries,
                    1-rightmargin-0.03,1-p1topmargin-0.03])
    # extra info box parameters
    if infoleft is None: infoleft = leftmargin+0.05
    if infotop is None: infotop = 1-p1topmargin-0.1
    # marker properties for data
    markerstyle = 20
    markercolor = 1
    markersize = 0.75

    ### style operations on theory histograms
    for i,hist in enumerate(theoryhists):
        if hist is not None:
            hist.SetLineWidth(2)
            hist.SetLineColor(colorlist[i])

    ### style operations on uncertainty histograms
    if systhists is not None:
        for i,hist in enumerate(systhists):
            if hist is not None:
                hist.SetLineWidth(0)
                hist.SetMarkerStyle(0)
                hist.SetFillStyle(3254)
                hist.SetFillColor(colorlist[i])

    ### style operations on data histogram
    datahist.SetMarkerStyle(markerstyle)
    datahist.SetMarkerColor(markercolor)
    datahist.SetMarkerSize(markersize)
    datahist.SetLineColor(markercolor)
    statdatahist.SetMarkerSize(0)
    statdatahist.SetLineColor(ROOT.kRed)

    ### make ratio histograms
    ratiohistlist = [h.Clone() for h in theoryhists]
    ratiodatahist = datahist.Clone()
    ratiostatdatahist = statdatahist.Clone()
    allratiohists = ratiohistlist+[ratiodatahist]+[ratiostatdatahist]
    ratiosysthists = None
    if systhists is not None:
        ratiosysthists = []
        for hist in systhists:
            if hist is not None: ratiosysthists.append(hist.Clone())
            else: ratiosysthists.append(None)
        allratiohists += ratiosysthists
    for hist in allratiohists:
        for j in range(0,hist.GetNbinsX()+2):
            scale = datahist.GetBinContent(j)
            if scale<1e-12:
                hist.SetBinContent(j,0)
                hist.SetBinError(j,10)
            else:
                hist.SetBinContent(j,hist.GetBinContent(j)/scale)
                hist.SetBinError(j,hist.GetBinError(j)/scale)
 
    ### make legend for upper plot and add all histograms
    legend = ROOT.TLegend(plegendbox[0],plegendbox[1],plegendbox[2],plegendbox[3])
    legend.SetNColumns(1)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(10*legendfont+3)
    legend.SetBorderSize(1)
    legend.AddEntry(datahist,'Data',"pe1")
    for i,hist in enumerate(theoryhists):
        label = hist.GetTitle()
        if labellist is not None: label = labellist[i]
        legend.AddEntry(hist,label,"l")

    ### make canvas and pads
    c1 = ROOT.TCanvas("c1","c1")
    c1.SetCanvasSize(cwidth,cheight)
    pad1 = ROOT.TPad("pad1","",0.,rfrac,1.,1.)
    pad1.SetTopMargin(p1topmargin)
    pad1.SetBottomMargin(p1bottommargin)
    pad1.SetLeftMargin(leftmargin)
    pad1.SetRightMargin(rightmargin)
    pad1.SetTicks(1,1)
    pad1.SetFrameLineWidth(2)
    pad1.SetGrid()
    pad1.Draw()
    pad2 = ROOT.TPad("pad2","",0.,0.,1.,rfrac)
    pad2.SetTopMargin(p2topmargin)
    pad2.SetBottomMargin(p2bottommargin)
    pad2.SetLeftMargin(leftmargin)
    pad2.SetRightMargin(rightmargin)
    pad2.SetTicks(1,1)
    pad2.SetFrameLineWidth(2)
    pad2.SetGrid()
    pad2.Draw()

    ### make upper part of the plot
    pad1.cd()

    # get x-limits (for later use)
    nbins = datahist.GetNbinsX()
    xlims = (datahist.GetBinLowEdge(1),
             datahist.GetBinLowEdge(nbins)+datahist.GetBinWidth(nbins))
    # get and set y-limits
    (totmin,totmax) = ht.getminmax(theoryhists+[datahist])
    # in case of log scale
    if logy:
        pad1.SetLogy()
        if ylims is None: ylims = (totmin*yminlogfactor, totmax*ymaxlogfactor)
    # in case of lin scale
    else:
        if ylims is None: ylims = (0.,totmax*ymaxlinfactor)
    if yminzero and ylims[0]<0: ylims = (0.,ylims[1])
    theoryhists[0].SetMaximum(ylims[1])
    theoryhists[0].SetMinimum(ylims[0])

    # X-axis layout
    xax = theoryhists[0].GetXaxis()
    xax.SetNdivisions(5,4,0,ROOT.kTRUE)
    xax.SetLabelSize(0)
    # Y-axis layout
    yax = theoryhists[0].GetYaxis()
    yax.SetMaxDigits(3)
    yax.SetNdivisions(8,4,0,ROOT.kTRUE)
    yax.SetLabelFont(10*labelfont+3)
    yax.SetLabelSize(labelsize)
    if yaxtitle is not None:
        yax.SetTitle(yaxtitle)
        yax.SetTitleFont(10*axtitlefont+3)
        yax.SetTitleSize(axtitlesize)
        yax.SetTitleOffset(ytitleoffset)

    # draw histograms
    theoryhists[0].Draw(drawoptions)
    if systhists is not None:
        for hist in systhists:
            if hist is not None: hist.Draw("same e2")
    for hist in theoryhists:
        hist.Draw("same "+drawoptions)
    datahist.Draw("pe e1 x0 same")
    statdatahist.Draw("pe e1 x0 same")
    if dolegend:
        legend.Draw("same")
    ROOT.gPad.RedrawAxis()

    # draw header
    pt.drawLumi(pad1, extratext=extracmstext, lumitext=lumitext, rfrac=rfrac)

    # draw extra info
    tinfo = ROOT.TLatex()
    tinfo.SetTextFont(10*infofont+3)
    tinfo.SetTextSize(infosize)
    for i,info in enumerate(extrainfos):
        vspace = 0.07*(float(infosize)/20)
        tinfo.DrawLatexNDC(infoleft,infotop-(i+1)*vspace, info)

    ### make the lower part of the plot
    pad2.cd()
    xax = ratiohistlist[0].GetXaxis()
    xax.SetNdivisions(5,4,0,ROOT.kTRUE)
    xax.SetLabelSize(labelsize)
    xax.SetLabelFont(10*labelfont+3)
    if xaxtitle is not None:
        xax.SetTitle(xaxtitle)
        xax.SetTitleFont(10*axtitlefont+3)
        xax.SetTitleSize(axtitlesize)
        xax.SetTitleOffset(xtitleoffset)
    # Y-axis layout
    yax = ratiohistlist[0].GetYaxis()
    if ratiorange==None: ratiorange = (0,1.999)
    yax.SetRangeUser(ratiorange[0],ratiorange[1]);
    yax.SetMaxDigits(3)
    yax.SetNdivisions(4,5,0,ROOT.kTRUE)
    yax.SetLabelFont(10*labelfont+3)
    yax.SetLabelSize(labelsize)
    yax.SetTitle('Ratio')
    yax.SetTitleFont(10*axtitlefont+3)
    yax.SetTitleSize(axtitlesize)
    yax.SetTitleOffset(ytitleoffset)

    # draw objects
    ratiohistlist[0].Draw(drawoptions)
    if ratiosysthists is not None:
        for hist in ratiosysthists:
            if hist is not None: hist.Draw("same e2")
    for hist in ratiohistlist:
        hist.Draw("same "+drawoptions)
    ratiodatahist.Draw("pe e1 x0 same")
    ratiostatdatahist.Draw("pe e1 x0 same")
    ROOT.gPad.RedrawAxis()

    # make and draw unit ratio line
    xmax = theoryhists[0].GetXaxis().GetBinUpEdge(theoryhists[0].GetNbinsX())
    xmin = theoryhists[0].GetXaxis().GetBinLowEdge(1)
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineStyle(2)
    line.Draw("same")

    c1.SaveAs(figname.replace('.png','')+'.png')
    #c1.SaveAs(figname.replace('.png','')+'.eps')
    #c1.SaveAs(figname.replace('.png','')+'.pdf')
