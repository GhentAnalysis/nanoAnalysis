#!/usr/bin/env python

#####################
# make impact plots #
#####################
# documentation: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts.
# this script was copied from CombineHarvester/CombineTools/scripts/plotImpacts.py
# with local modifications to allow multi-POI plotting (not supported in the builtin version).

import ROOT
import math
import json
import argparse
import CombineHarvester.CombineTools.plotting as plot
import CombineHarvester.CombineTools.combine.rounding as rounding

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='input json file')
parser.add_argument('--output', '-o', help='name of the output file to create')
parser.add_argument('--translate', '-t', help='JSON file for remapping of parameter names')
parser.add_argument('--units', default=None, help='Add units to the best-fit parameter value')
parser.add_argument('--per-page', type=int, default=30, help='Number of parameters to show per page')
parser.add_argument('--max-pages', type=int, default=None, help='Maximum number of pages to write')
parser.add_argument('--height', type=int, default=600, help='Canvas height, in pixels')
parser.add_argument('--left-margin', type=float, default=0.4, help='Left margin, expressed as a fraction')
parser.add_argument('--label-size', type=float, default=0.021, help='Parameter name label size')
parser.add_argument('--cms-label', default='Internal', help='Label next to the CMS logo')
parser.add_argument('--transparent', action='store_true', help='Draw areas as hatched lines instead of solid')
parser.add_argument('--checkboxes', action='store_true', help='Draw an extra panel with filled checkboxes')
parser.add_argument('--blind', action='store_true', help='Do not print best fit signal strength')
parser.add_argument('--color-groups', default=None, help='Comma separated list of GROUP=COLOR')
parser.add_argument("--pullDef",  default=None, help="Choose the definition of the pull, see HiggsAnalysis/CombinedLimit/python/calculate_pulls.py for options")
parser.add_argument('--POIs', default=None, nargs='+', help='Specify the POIs to draw. Default: first POI in input file.')
args = parser.parse_args()

if args.transparent:
    print 'plotImpacts.py: --transparent is now always enabled, the option will be removed in a future update'

if args.checkboxes:
    args.checkboxes = False
    msg = 'WARNING: option checkboxes overridden to False.'
    print(msg)

externalPullDef = False
if args.pullDef is not None:
    externalPullDef = True
    import HiggsAnalysis.CombinedLimit.calculate_pulls as CP


def Translate(name, ndict):
    return ndict[name] if name in ndict else name


def GetRounded(nom, e_hi, e_lo):
    if e_hi < 0.0:
        e_hi = 0.0
    if e_lo < 0.0:
        e_lo = 0.0
    # original code:
    rounded = rounding.PDGRoundAsym(nom, e_hi if e_hi != 0.0 else 1.0, e_lo if e_lo != 0.0 else 1.0)
    # local change: increase significant digits
    #rounded = (str(int(round(nom * 1000))),
    #                [[str(int(round(e_hi*1000))),str(int(round(e_lo*1000)))]],
    #                -3 )
    s_nom = rounding.downgradePrec(rounded[0],rounded[2])
    s_hi = rounding.downgradePrec(rounded[1][0][0],rounded[2]) if e_hi != 0.0 else '0'
    s_lo = rounding.downgradePrec(rounded[1][0][1],rounded[2]) if e_lo != 0.0 else '0'
    return (s_nom, s_hi, s_lo)

# Dictionary to translate parameter names
translate = {}
if args.translate is not None:
    with open(args.translate) as jsonfile:
        translate = json.load(jsonfile)

# Load the json output of combineTool.py -M Impacts
data = {}
with open(args.input) as jsonfile:
    data = json.load(jsonfile)

# Set the global plotting style
plot.ModTDRStyle(l=args.left_margin, b=0.10, width=(900 if args.checkboxes else 700), height=args.height)

# Find POIs to plot
allPOIs = [ele['name'] for ele in data['POIs']]
POIs = [allPOIs[0]]
if args.POIs is None:
    msg = 'The --POI argument was not specified;'
    msg += ' making impact plot for first POI found in input file,'
    msg += ' named {}'.format(POIs[0])
else:
    POIs = []
    for POI in args.POIs:
        if not POI in allPOIs:
            msg = 'ERROR: requested POI {} not found in input file.'.format(POI)
            raise Exception(msg)
        POIs.append(POI)
POIs_info = {}
POIs_fit = {}
for ele in data['POIs']:
    if ele['name'] in POIs:
        POIs_info[ele['name']] = ele
        POIs_fit[ele['name']] = ele['fit']
firstPOI = POIs[0]

# Sort parameters by largest summed absolute impact on all POIs
data['params'].sort(key=lambda x: sum( [abs(x['impact_%s' % POI]) for POI in POIs] ) , reverse=True)

if args.checkboxes:
    cboxes = data['checkboxes']

# Set the number of parameters per page (show) and the number of pages (n)
show = args.per_page
n = int(math.ceil(float(len(data['params'])) / float(show)))
if args.max_pages is not None and args.max_pages > 0:
    n = args.max_pages

colors = {
    'Gaussian': 1,
    'Poisson': 8,
    'AsymmetricGaussian': 9,
    'Unconstrained': 39,
    'Unrecognised': 2
}
color_hists = {}
color_group_hists = {}

if args.color_groups is not None:
    color_groups = {
        x.split('=')[0]: int(x.split('=')[1]) for x in args.color_groups.split(',')
    }

seen_types = set()

for name, col in colors.iteritems():
    color_hists[name] = ROOT.TH1F()
    plot.Set(color_hists[name], FillColor=col, Title=name)

if args.color_groups is not None:
    for name, col in color_groups.iteritems():
        color_group_hists[name] = ROOT.TH1F()
        plot.Set(color_group_hists[name], FillColor=col, Title=name)

for page in xrange(n):
    canv = ROOT.TCanvas(args.output, args.output)
    n_params = len(data['params'][show * page:show * (page + 1)])
    pdata = data['params'][show * page:show * (page + 1)]
    print '>> Doing page %i, have %i parameters' % (page, n_params)

    # make the horizontal white/gray bands with parameter ranking number
    boxes = []
    for i in xrange(n_params):
        y1 = ROOT.gStyle.GetPadBottomMargin()
        y2 = 1. - ROOT.gStyle.GetPadTopMargin()
        h = (y2 - y1) / float(n_params)
        y1 = y1 + float(i) * h
        y2 = y1 + h
        box = ROOT.TPaveText(0, y1, 1, y2, 'NDC')
        plot.Set(box, TextSize=0.02, BorderSize=0, FillColor=0, TextAlign=12, Margin=0.005)
        if i % 2 == 0:
            box.SetFillColor(18)
        box.AddText('%i' % (n_params - i + page * show))
        box.Draw()
        boxes.append(box)

    # Crate and style the pads
    if args.checkboxes:
        pads = plot.MultiRatioSplitColumns([0.54, 0.24], [0., 0.], [0., 0.])
        pads[2].SetGrid(1, 0)
    else:
        margin = 0.015
        splitpoints = [0.7]
        if len(POIs)>1:
            splitpoints = [0.4]
            for i in range(1, len(POIs)):
                width = 0.6 / len(POIs) 
                width -= margin
                splitpoints.append(width)
        pads = plot.MultiRatioSplitColumns(splitpoints, [margin]*len(splitpoints), [margin]*len(splitpoints))
    pads[0].SetGrid(1, 0)
    pads[0].SetTickx(1)
    pads[0].SetLeftMargin(0.2)
    for i in range(1, len(pads)):
        pads[i].SetGrid(1, 0)
        pads[i].SetTickx(1)

    h_pulls = ROOT.TH2F("pulls", "pulls", 6, -2.9, 2.9, n_params, 0, n_params)
    g_pulls = ROOT.TGraphAsymmErrors(n_params)
    g_impacts_hi = {}
    for POI in POIs: g_impacts_hi[POI] = ROOT.TGraphAsymmErrors(n_params)
    g_impacts_lo = {}
    for POI in POIs: g_impacts_lo[POI] = ROOT.TGraphAsymmErrors(n_params)
    g_check = ROOT.TGraphAsymmErrors()
    g_check_i = 0

    max_impact = 0.

    text_entries = [] # only for unconstrained parameters
    redo_boxes = [] # only for unconstrained parameters

    # loop over all parameters
    for p in xrange(n_params):
        i = n_params - (p + 1)
        pre = pdata[p]['prefit']
        fit = pdata[p]['fit']
        tp = pdata[p]['type']
        seen_types.add(tp)
        # case of all types except 'Unconstrained'
        if pdata[p]['type'] != 'Unconstrained':
            pre_err_hi = (pre[2] - pre[1])
            pre_err_lo = (pre[1] - pre[0])

            if externalPullDef:
                fit_err_hi = (fit[2] - fit[1])
                fit_err_lo = (fit[1] - fit[0])
                pull, pull_hi, pull_lo = CP.returnPullAsym(args.pullDef,fit[1],pre[1],fit_err_hi,pre_err_hi,fit_err_lo,pre_err_lo)
            else:
                pull = fit[1] - pre[1]
                pull = (pull/pre_err_hi) if pull >= 0 else (pull/pre_err_lo)
                pull_hi = fit[2] - pre[1]
                pull_hi = (pull_hi/pre_err_hi) if pull_hi >= 0 else (pull_hi/pre_err_lo)
                pull_hi = pull_hi - pull
                pull_lo = fit[0] - pre[1]
                pull_lo = (pull_lo/pre_err_hi) if pull_lo >= 0 else (pull_lo/pre_err_lo)
                pull_lo =  pull - pull_lo

            g_pulls.SetPoint(i, pull, float(i) + 0.5)
            g_pulls.SetPointError(i, pull_lo, pull_hi, 0., 0.)
        else:
            # special case for 'Unconstrained' parameter:
            # hide the pull point and make special text entry
            g_pulls.SetPoint(i, 0., 9999.)
            y1 = ROOT.gStyle.GetPadBottomMargin()
            y2 = 1. - ROOT.gStyle.GetPadTopMargin()
            x1 = ROOT.gStyle.GetPadLeftMargin()
            h = (y2 - y1) / float(n_params)
            y1 = y1 + ((float(i)+0.5) * h)
            x1 = x1 + (1 - pads[0].GetRightMargin() -x1)/2.
            s_nom, s_hi, s_lo = GetRounded(fit[1], fit[2] - fit[1], fit[1] - fit[0])
            text_entries.append((x1, y1, '%s^{#plus%s}_{#minus%s}' % (s_nom, s_hi, s_lo)))
            redo_boxes.append(i)
        if args.checkboxes:
            pboxes = pdata[p]['checkboxes']
            for pbox in pboxes:
                cboxes.index(pbox)
                g_check.SetPoint(g_check_i, cboxes.index(pbox) + 0.5, float(i) + 0.5)
                g_check_i += 1
        # fill impact graphs
        for POI in POIs:
            g_impacts_hi[POI].SetPoint(i, 0, float(i) + 0.5)
            g_impacts_lo[POI].SetPoint(i, 0, float(i) + 0.5)
            imp = pdata[p][POI]
            g_impacts_hi[POI].SetPointError(i, 0, imp[2] - imp[1], 0.5, 0.5)
            g_impacts_lo[POI].SetPointError(i, imp[1] - imp[0], 0, 0.5, 0.5)
            max_impact = max(
                max_impact, abs(imp[1] - imp[0]), abs(imp[2] - imp[1]))
        # find correct color for parameter name
        col = colors.get(tp, 2)
        if args.color_groups is not None and len(pdata[p]['groups']) >= 1:
            for p_grp in pdata[p]['groups']:
                if p_grp in color_groups:
                    col = color_groups[p_grp]
                    break
        # set parameter name as bin label in pull histogram
        h_pulls.GetYaxis().SetBinLabel(
            i + 1, ('#color[%i]{%s}'% (col, Translate(pdata[p]['name'], translate))))

    # Style and draw the pulls histo
    if externalPullDef:
        plot.Set(h_pulls.GetXaxis(), TitleSize=0.04, LabelSize=0.03, Title=CP.returnTitle(args.pullDef))
    else:
        plot.Set(h_pulls.GetXaxis(), TitleSize=0.03, LabelSize=0.02, Title='(#hat{#theta}-#theta_{0})/#Delta#theta')
    plot.Set(h_pulls.GetYaxis(), LabelSize=args.label_size, TickLength=0.0)
    h_pulls.GetYaxis().LabelsOption('v')
    h_pulls.Draw()

    # write info for unconstrained parameters
    for i in redo_boxes:
        newbox = boxes[i].Clone()
        newbox.Clear()
        newbox.SetY1(newbox.GetY1()+0.005)
        newbox.SetY2(newbox.GetY2()-0.005)
        newbox.SetX1(ROOT.gStyle.GetPadLeftMargin()+0.001)
        newbox.SetX2(0.7-0.001)
        newbox.Draw()
        boxes.append(newbox)
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.02)
    latex.SetTextAlign(22)
    for entry in text_entries:
        latex.DrawLatex(*entry)

    # Go to the other pad and draw the impacts histos
    h_impacts = {}
    for pad, POI in zip(pads[1:], POIs):
        pad.cd()
        if max_impact == 0.: max_impact = 1E-6  # otherwise the plotting gets screwed up
        impacthist = ROOT.TH2F(
            "impacts", "impacts", 6, -max_impact * 1.1, max_impact * 1.1, n_params, 0, n_params)
        plot.Set(impacthist.GetXaxis(), LabelSize=0.02, TitleSize=0.03, Ndivisions=505, Title=
            '#Delta#hat{%s}' % (Translate(POI, translate)))
        plot.Set(impacthist.GetYaxis(), LabelSize=0, TickLength=0.0)
        impacthist.Draw()
        h_impacts[POI] = impacthist

        if args.checkboxes:
            pads[2].cd()
            h_checkboxes = ROOT.TH2F(
                "checkboxes", "checkboxes", len(cboxes), 0, len(cboxes), n_params, 0, n_params)
            for i, cbox in enumerate(cboxes):
                h_checkboxes.GetXaxis().SetBinLabel(i+1, Translate(cbox, translate))
            plot.Set(h_checkboxes.GetXaxis(), LabelSize=0.03, LabelOffset=0.002)
            h_checkboxes.GetXaxis().LabelsOption('v')
            plot.Set(h_checkboxes.GetYaxis(), LabelSize=0, TickLength=0.0)
            h_checkboxes.Draw()
            # g_check.SetFillColor(ROOT.kGreen)
            g_check.Draw('PSAME')

    # Back to the first pad to draw the pulls graph
    pads[0].cd()
    plot.Set(g_pulls, MarkerSize=0.8, LineWidth=2)
    g_pulls.Draw('PSAME')

    # And back to the other pads to draw the impacts graphs
    alpha = 0.7
    lo_color = {
        'default': 38,
        'hesse': ROOT.kOrange - 3,
        'robust': ROOT.kGreen + 1
    }
    hi_color = {
        'default': 46,
        'hesse': ROOT.kBlue,
        'robust': ROOT.kAzure - 5
    }
    method = 'default'
    if 'method' in data and data['method'] in lo_color:
        method = data['method']
    for pad, POI in zip(pads[1:], POIs):
        pad.cd()
        g_impacts_hi[POI].SetFillColor(plot.CreateTransparentColor(hi_color[method], alpha))
        g_impacts_hi[POI].Draw('2SAME')
        g_impacts_lo[POI].SetFillColor(plot.CreateTransparentColor(lo_color[method], alpha))
        g_impacts_lo[POI].Draw('2SAME')
        pad.RedrawAxis()

    # make a legend
    legend = ROOT.TLegend(0.02, 0.02, 0.20, 0.08, '', 'NBNDC')
    legend.SetNColumns(1)
    legend.AddEntry(g_pulls, 'Pull', 'LP')
    legend.AddEntry(g_impacts_hi[firstPOI], '+1#sigma Impact', 'F')
    legend.AddEntry(g_impacts_lo[firstPOI], '-1#sigma Impact', 'F')
    legend.Draw()

    # make another legend if there are different types of parameters
    leg_width = pads[0].GetLeftMargin() - 0.01
    if args.color_groups is not None:
        legend2 = ROOT.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
        legend2.SetNColumns(2)
        for name, h in color_group_hists.iteritems():
            legend2.AddEntry(h, Translate(name, translate), 'F')
        legend2.Draw()
    elif len(seen_types) > 1:
        legend2 = ROOT.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
        legend2.SetNColumns(2)
        for name, h in color_hists.iteritems():
            if name == 'Unrecognised': continue
            legend2.AddEntry(h, name, 'F')
        legend2.Draw()

    # draw CMS header and pad titles
    plot.DrawCMSLogo(pads[0], 'CMS', args.cms_label, 0, 0., 0., 0., cmsTextSize=0.5)
    for pad, POI in zip(pads[1:], POIs):
        fit = POIs_fit[POI]
        s_nom, s_hi, s_lo = GetRounded(fit[1], fit[2] - fit[1], fit[1] - fit[0])
        if not args.blind:
            plot.DrawTitle(pad, '#hat{%s} = %s^{#plus%s}_{#minus%s}%s' % (
                Translate(POI, translate), s_nom, s_hi, s_lo,
                    '' if args.units is None else ' '+args.units), 3, 0.27, textSize=0.4)
    extra = ''
    if page == 0:
        extra = '('
    if page == n - 1:
        extra = ')'
    canv.Print('.pdf%s' % extra)
