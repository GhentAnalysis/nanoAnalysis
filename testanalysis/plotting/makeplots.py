###############################
# Make basic MC vs data plots #
###############################
# Note: this script does not take into account systematic uncertainties,
#       only nominal histograms and statistical uncertainties.
# Input histograms:
#   The input histograms are supposed to be contained in a single ROOT file.
#   The naming of the histograms should be <process name>_<region>_<variable name>.
#   Optionally, their names can be <process name>_<region>_<variable name>_<systematic>,
#   but in that case only the nominal ones are used and others are ignored.

# import python modules
import sys
import os
import argparse
#from pathlib import Path
# import framework modules
#sys.path.append(str(Path(__file__).parents[2]))
#import tools.histtools as ht
#from tools.variabletools import read_variables
#import plotting.histplotter as hp
#from constants.luminosities import lumidict
sys.path.append(os.path.abspath('../../tools'))
import histtools as ht
from variabletools import read_variables
sys.path.append(os.path.abspath('../../plotting'))
import histplotter as hp
sys.path.append(os.path.abspath('../../constants'))
from luminosities import lumidict
# import local modules
import colors
import infodicts


if __name__=="__main__":

  # parse arguments
  parser = argparse.ArgumentParser(description='Make plots')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-y', '--year', required=True)
  parser.add_argument('-r', '--region', required=True)
  parser.add_argument('-v', '--variables', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
  parser.add_argument('--colormap', default='default')
  parser.add_argument('--signals', default=None, nargs='+')
  parser.add_argument('--extracmstext', default='Preliminary')
  parser.add_argument('--unblind', action='store_true')
  parser.add_argument('--dolog', action='store_true')
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))
    
  # read all histogram names
  histnames = ht.loadallhistnames(args.inputfile)

  # check if some of the histogram names end with '_nominal',
  # and if so, select only those histogram names
  hassystematics = False
  for histname in histnames:
    if histname.endswith('_nominal'):
      hassystematics = True
      break
  if hassystematics:
    histnames = [histname for histname in histnames if histname.endswith('_nominal')]

  # load (nominal) histograms 
  histlist = ht.loadhistogramlist(args.inputfile, histnames)
  ht.cliphistograms(histlist)

  # re-set histogram name and title
  # (could diverge from key name if only key name was modified for speed)
  for hist, histname in zip(histlist, histnames):
    hist.SetName(histname)
    hist.SetTitle(histname.split('_')[0])

  # select histograms
  histlist = ht.selecthistograms(histlist,mustcontainall=[args.region])[1]
  if len(histlist)==0:
    raise Exception('ERROR: histogram list is empty, cannot make plots.')

  # make output directory
  if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

  # get a printable version of the region name
  regiondict = infodicts.get_region_dict()
  if args.region in regiondict.keys():
    regionname = regiondict[args.region]
  else:
    print('WARNING: region {} not found in region dict,'.format(args.region),
          ' will write raw region name on plot.')
    regionname = args.region

  # read variables
  variables = read_variables( args.variables )

  # loop over variables
  for var in variables:
    varname = var.name
    axtitle = var.axtitle
    unit = var.unit
    
    # select histograms
    thishists = ht.selecthistograms(histlist, mustcontainall=['_'+varname])[1]
    # additional selections for overlapping histogram names
    thishists = ([hist for hist in thishists if 
                  (hist.GetName().endswith(varname) or varname+'_' in hist.GetName())])
    if len(thishists)==0:
      print('ERROR: histogram list for variable {} is empty,'.format(varname)
            +' skipping this variable.')
      continue

    # printouts for testing
    #for hist in thishists: print(hist.GetName())
    
    # find data and sim histograms
    datahists = []
    simhists = []
    for hist in thishists:
      if( hist.GetName().startswith('data') or hist.GetName().startswith('Data') ): 
        datahists.append(hist)
      else: simhists.append(hist)
    if len(datahists)==0:
      print('WARNING: no data histogram found, plotting simulation only.')
      datahist = simhists[0].Clone()
      args.unblind = False
    elif len(datahists)!=1:
      msg = 'ERROR: expecting one data histogram'
      msg += ' but found {}:\n'.format(len(datahists))
      for datahist in datahists: msg += '  - {}\n'.format(datahist.GetName())
      msg.strip('\n')
      raise Exception(msg)
    else: datahist = datahists[0]

    # blind data histogram
    if not args.unblind:
      for i in range(0,datahist.GetNbinsX()+2):
        datahist.SetBinContent(i, 0)
        datahist.SetBinError(i, 0)

    # set plot properties
    xaxtitle = axtitle
    if( axtitle is not None and unit is not None ):
      xaxtitle += ' ({})'.format(unit)
    yaxtitle = 'Events'
    outfile = os.path.join(args.outputdir, varname)
    if not args.year in lumidict.keys():
      print('WARNING: year {} not recognized,'.format(args.year)
            +' will not write lumi header.')
    lumi = lumidict.get(args.year,None)
    colormap = colors.getcolormap(style=args.colormap)
    extrainfos = []
    extrainfos.append( args.year )
    extrainfos.append( regionname )
    xlabels = None
    labelsize = None
    if( var.iscategorical and var.xlabels is not None ):
        xlabels = var.xlabels
        labelsize = 15

    # make the plot
    hp.plotdatavsmc(outfile, datahist, simhists,
	    xaxtitle=xaxtitle,
	    yaxtitle='Number of events',
	    colormap=colormap,
	    lumi=lumi, extracmstext=args.extracmstext,
            extrainfos=extrainfos, infosize=15,
            binlabels=xlabels, labelsize=labelsize,
            signals=args.signals )

    if args.dolog:
      # make plot in log scale
      outfile = os.path.join(args.outputdir, varname)+'_log'
      hp.plotdatavsmc(outfile, datahist, simhists,
            xaxtitle=xaxtitle,
            yaxtitle='Number of events',
            yaxlog=True,
            colormap=colormap,
            lumi=lumi, extracmstext=args.extracmstext,
            extrainfos=extrainfos, infosize=15,
            binlabels=xlabels, labelsize=labelsize,
            signals=args.signals )
