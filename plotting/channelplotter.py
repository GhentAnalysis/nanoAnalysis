############################################################
# Make a plot of results for different channels to compare #
############################################################

import sys
import os
import math
import ROOT
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
import plotting.plottools as pt

def readchanneltxt(path_to_txtfile):
    ### read a txt file and return a list of channels that can be used for plotchannels
    # expected formatting of txt file: name label central -statdown +statup -totdown +totup
    # note: now more flexible allowing for any number of entries per channel,
    #       as long as this number is the same for each channel.
    # note: label may not contain spaces, but underscores can be used instead
    # note: lines in the file starting with # will be ignored
    channels = []
    nentries = -1
    f = open(path_to_txtfile,'r')
    for line in f:
        if line[0]=='#': continue
        line = line.strip(' ')
        channel = line.split(' ')
        if nentries>0 and len(channel)!=nentries:
            raise Exception('ERROR: inconsistent number of entries per channel'
                            +' while reading input file {}'.format(path_to_txtfle))
        nentries = len(channel)
        # replace underscores by spaces except when a subscript is implied
        channel[1] = channel[1].replace('_{','blablatemp')
        channel[1] = channel[1].replace('_',' ')
        channel[1] = channel[1].replace('blablatemp','_{')
        for i in range(len(channel))[2:]: channel[i] = float(channel[i].strip('+-'))
        channels.append(channel)
    return channels

def get_xaxrange(channels):
    ### automatically determine a suitable x-axis range
    ups = []
    downs = []
    for cnum,channel in enumerate(channels):
        central = float(channel[2])
        totdown = float(channel[5])
        totup = float(channel[6])
        ups.append(central+totup)
        downs.append(central-totdown)
    cmin = min(downs)
    cmax = max(ups)
    margin = (cmax-cmin)/5.
    return (cmin-5*margin, cmax+margin)

def plotchannels(channels, outfile,
                showvalues=True, font=None, fontsize=None,
                xaxtitle=None, yaxtitle=None, lumi=None,
                xaxcentral=None, xaxlinecoords=[], yaxlinecoords=[],
                xaxrange=None, legendbox=None,
                extracmstext=''):
    ### make the plot of different channels
    ### arguments:
    # - channels: list of lists of the form [name,label,central,statdown,statup,totdown,totup]
    #             (note: name attribute is not used here)
    # - outfile: name of the output file to make (should not contain the file type extension!)
    # - showvalues: boolean whether to plot the numerical values or only the labels
    # - xaxtitle, yaxtitle are axis titles (default: no title)
    # - lumi is luminosity value to put in header
    # - xaxcentral is x-axis value where to draw a vertical line (default: no vertical line)
    # - xaxlinecoords: list of x-axis coordinates where to draw additional vertical lines
    # - xaxrange is range of x axis (default: 0-1.5)

    pt.setTDRstyle()
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)

    # parse arguments
    if font is None: font = 43 
    if fontsize is None: fontsize = 20
    lumi = str(lumi)
    dolumi = True
    if lumi=='': dolumi=False
    if xaxrange is None: xaxrange=get_xaxrange(channels)
    docentral = True
    if xaxcentral is None: 
        docentral = False
        xaxcentral = 0
    legendnrows = 1.5
    if legendbox is None: 
        legendbox = [0.15, 0.87, 0.9, 0.92]

    # intializations
    c1 = ROOT.TCanvas("c1")
    pad1 = ROOT.TPad("pad1","",0.,0.,1.,1.)
    pad1.SetBottomMargin(0.15)
    pad1.SetLeftMargin(0.1)
    pad1.SetTopMargin(0.06)
    #pad1.SetGridx(1)
    #pad1.SetGridy(1)
    pad1.Draw()
    pad1.cd()
    nchannels = len(channels)
    points = ROOT.TGraphAsymmErrors(nchannels)
    points_statonly = ROOT.TGraphAsymmErrors(nchannels)
    labels = []
    mg = ROOT.TMultiGraph()
    mg.Add(points)
    mg.Add(points_statonly)   

    # loop over channels and make points and appropriate labels
    for cnum,channel in enumerate(channels):
        label = str(channel[1])
        central = float(channel[2])
        statdown = float(channel[3])
        statup = float(channel[4])
        totdown = float(channel[5])
        totup = float(channel[6])
        if totup>statup:
            systup = math.sqrt(totup**2-statup**2)
        else:
            msg = 'WARNING in channelplotter: stat-only up error > total up error'
            msg += ' for channel {} ({})'.format(channel[0],label)
            print(msg)
            systup = totup
            statup = 0
        if totdown>statdown:
            systdown = math.sqrt(totdown**2-statdown**2)
        else:
            msg = 'WARNING in channelplotter: stat-only down error > total down error'
            msg += ' for channel {} ({})'.format(channel[0],label)
            print(msg)
            systdown = totdown
            statdown = 0

        points.SetPoint(cnum, central, cnum+1)
        points.SetPointError(cnum, totdown, totup, 0.5, 0.5)
        
        points_statonly.SetPoint(cnum, central, cnum+1)
        points_statonly.SetPointError(cnum, statdown, statup, 0.5, 0.5)

        if showvalues:
            centraltxt = '{:.2f}'.format(central)
            statuptxt = '{:.2f}'.format(statup)
            statdowntxt = '{:.2f}'.format(statdown)
            systuptxt = '{:.2f}'.format(systup)
            systdowntxt = '{:.2f}'.format(systdown)
            valtxt = centraltxt
            valtxt += '^{+'+statuptxt+'}_{-'+statdowntxt+'} (stat)'
            valtxt += '^{+'+systuptxt+'}_{-'+systdowntxt+'} (syst)'
            label += ': '+valtxt

        labels.append(label)

    # set style for points
    points.SetFillStyle(1001)
    points.SetFillColor(ROOT.kOrange-4)
    points.SetMarkerStyle(20)

    # set style for statonly points
    points_statonly.SetFillStyle(1001)
    points_statonly.SetFillColor(ROOT.kCyan-4)
    points_statonly.SetMarkerStyle(20)
    
    # draw first time to have well-defined axes
    mg.Draw("A2")
    pad1.Update()
    c1.Update()

    # x-axis layout
    mg.GetXaxis().SetLimits(xaxrange[0],xaxrange[1])
    if xaxtitle is not None:
        mg.GetXaxis().SetTitle(xaxtitle)
        mg.GetXaxis().SetTitleSize(0.05)
    mg.GetXaxis().SetLabelSize(0.04)
    
    # y-axis layout
    yax = mg.GetYaxis()
    if yaxtitle is not None:
        yax.SetTitle(yaxtitle)
        yax.SetTitleSize(0.05)
        yax.SetTitleOffset(0.5)
    yax.SetLabelSize(0)
    yax.SetTickLength(0)
    mg.SetMaximum(len(channels)+legendnrows)
    
    # unit line
    sm = ROOT.TLine(xaxcentral,0.5,xaxcentral,nchannels+0.5)
    # style 1: gray band
    #sm.SetLineColor(ROOT.kGray)
    #sm.SetLineWidth(3)
    #sm.SetLineStyle(1)
    # style 2: black dashes
    sm.SetLineColor(ROOT.kBlack)
    sm.SetLineWidth(1)
    sm.SetLineStyle(9)

    # additional lines
    vertlines = []
    for xcoord in xaxlinecoords:
        vertlines.append( ROOT.TLine(xcoord,0.3,xcoord,nchannels+0.7) )
        vertlines[-1].SetLineWidth(1)
        vertlines[-1].SetLineStyle(3)
    horlines = []
    for ycoord in yaxlinecoords:
        horlines.append( ROOT.TLine(xaxrange[0],ycoord,xaxrange[1],ycoord) )
        horlines[-1].SetLineWidth(1)
        horlines[-1].SetLineStyle(3)

    # legend
    leg = ROOT.TLegend(legendbox[0],legendbox[1],legendbox[2],legendbox[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(points,'Measurement',"p")
    leg.AddEntry(points_statonly, 'Stat.', "f")
    leg.AddEntry(points, 'Total', "f")
    leg.SetNColumns(3)

    # draw objects
    mg.Draw("P A2")
    if docentral: sm.Draw()
    for vertline in vertlines: vertline.Draw()
    for horline in horlines: horline.Draw()
    pad1.Update()
    c1.Update()
    for cnum,label in enumerate(labels):
        tex = ROOT.TLatex()
        tex.SetTextFont(font)
        tex.SetTextSize(fontsize)
        ypos = points.GetY()[cnum]
        yposndc = (ypos - ROOT.gPad.GetY1())/(ROOT.gPad.GetY2()-ROOT.gPad.GetY1())
        tex.DrawLatexNDC(0.15,yposndc,label)
    leg.Draw()
    lumitext = None if lumi is None else str(lumi)+' fb^{-1} (13 TeV)'
    pt.drawLumi(pad1, extratext=extracmstext, lumitext=lumitext, cms_in_grid=False,
                 cmstext_size_factor=0.6, cmstext_offset=0.01)

    pad1.Update()
    c1.Update()

    ### save the plot
    c1.SaveAs(outfile.replace('.png','')+'.png')
    c1.SaveAs(outfile.replace('.png','')+'.eps')
    c1.SaveAs(outfile.replace('.png','')+'.pdf')

if __name__=='__main__':

    if len(sys.argv)!=3:
        print('### ERROR ###: incorrect number of command line arguments')
        print('               need <inputfile> <outputfile>')
        sys.exit()

    txtfile = os.path.abspath(sys.argv[1])
    figfile = sys.argv[2]
    if not os.path.exists(txtfile):
        print('### ERROR ###: input file does not seem to exist...') 
        sys.exit()
    
    channels = readchanneltxt(txtfile)

    plotchannels(channels, figfile,
                 showvalues=True,
                 xaxcentral=1,
                 xaxtitle='Observed / Predicted')
    print('done!')
