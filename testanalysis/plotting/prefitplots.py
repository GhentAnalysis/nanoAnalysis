##############
# make plots #
##############

# The input histograms are supposed to be contained in a single root file.
# The naming of the histograms should be <process name>_<region>_<variable name>_<systematic>
# where the systematic is either "nominal" or a systematic name followed by "Up" or "Down".

# import python modules
import sys
import os
import ROOT
import argparse
# import framework modules
sys.path.append(os.path.abspath('../../tools'))
import histtools as ht
import listtools as lt
from variabletools import HistogramVariable
from variabletools import DoubleHistogramVariable
from variabletools import read_variables
from processinfo import ProcessInfoCollection, ProcessCollection
sys.path.append(os.path.abspath('../../plotting'))
import histplotter as hp
sys.path.append(os.path.abspath('../combine'))
from uncertaintytools import remove_systematics_default
from uncertaintytools import add_systematics_default
from uncertaintytools import remove_systematics_all
from uncertaintytools import add_systematics_dummy
# import local modules
import colors
import infodicts
sys.path.append(os.path.abspath('../tools'))
from histogramselection import select_histnames


if __name__=="__main__":

  # parse arguments
  parser = argparse.ArgumentParser(description='Make prefit plots')
  parser.add_argument('-i', '--inputfile', required=True, type=os.path.abspath)
  parser.add_argument('-y', '--year', required=True)
  parser.add_argument('-r', '--region', required=True)
  parser.add_argument('-p', '--processes', required=True, nargs='+',
                      help='List of process tags to take into account;'
                          +' use "all" to use all processes in the input file.')
  parser.add_argument('-v', '--variables', required=True, type=os.path.abspath,
                      help='Path to json file holding variable definitions.')
  parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
  parser.add_argument('-d', '--datatag', default='Data')
  parser.add_argument('-c', '--colormap', default='default')
  parser.add_argument('-s', '--signals', default=None, nargs='+')
  parser.add_argument('--includetags', default=None, nargs='+',
                      help='List of systematic tags to include')
  parser.add_argument('--excludetags', default=None, nargs='+',
                      help='List of systematic tags to exclude')
  parser.add_argument('--tags', default=None, nargs='+',
                      help='List of additional info to display on plot'
                          +' (e.g. simulation year or selection region).'
                          +' Use underscores for spaces.')
  parser.add_argument('--extracmstext', default='Preliminary')
  parser.add_argument('--splitvariable',default=None)
  parser.add_argument('--splitprocess',default=None)
  parser.add_argument('--unblind', default=False, action='store_true')
  parser.add_argument('--dolog', default=False, action='store_true')
  parser.add_argument('--rawsystematics', default=False, action='store_true',
                      help='Take the systematics from the input file without modifications'
                          +' (i.e. no disablings and no adding of norm uncertainties).')
  parser.add_argument('--dummysystematics', default=False, action='store_true',
                      help='Use dummy systematics (see uncertaintytools for details).')
  args = parser.parse_args()

  # print arguments
  print('Running with following configuration:')
  for arg in vars(args):
    print('  - {}: {}'.format(arg,getattr(args,arg)))
    
  # parse input file
  if not os.path.exists(args.inputfile):
    raise Exception('ERROR: requested to run on '+args.inputfile
                    +' but it does not seem to exist...')

  # parse the string with process tags
  doallprocesses = (len(args.processes)==1 and args.processes[0]=='all')

  # parse the variables
  varlist = read_variables(args.variables)
  variablenames = [v.name for v in varlist]

  # parse tags
  extratags = []
  if args.tags is not None: extratags = arg.tags
  extratags = [t.replace('_',' ') for t in extratags]

  # make the output directory
  if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

  # get a printable version of the region name
  regiondict = infodicts.get_region_dict()
  if args.region in regiondict.keys():
    regionname = regiondict[args.region]
  else:
    print('WARNING: region {} not found in region dict,'.format(args.region)
          +' will write raw region name on plot.')
    regionname = args.region

  # get a dictionary to match histogram titles to legend entries
  processdict = infodicts.get_process_dict()

  # get all relevant histograms
  print('Loading histogram names from input file...')
  histnames = select_histnames(args.inputfile,
    processes=args.processes,
    regions=[args.region],
    variablenames=variablenames,
    includesystematics=args.includetags,
    excludesystematics=args.excludetags,
    verbose=True)

  # print all histogram names (only for testing)
  #print('Found following histograms')
  #for histname in histnames:
  #  print('  - {}'.format(histname))

  # make a ProcessInfoCollection to extract information
  # (use first variable, assume list of processes, systematics etc.
  #  is the same for all variables)
  splittag = args.region+'_'+variablenames[0]
  print('Constructing ProcessInfoCollection using split tag "{}"'.format(splittag))
  PIC = ProcessInfoCollection.fromhistlist( histnames, splittag, datatag=args.datatag )

  # manage systematics (not yet needed here, but useful for printing the correct info)
  if( not args.rawsystematics and not args.dummysystematics ):
    _ = remove_systematics_default( PIC, year=args.year )
    _ = add_systematics_default( PIC, year=args.year )
  if args.dummysystematics:
    _ = remove_systematics_all( PIC )
    _ = add_systematics_dummy( PIC )
  #print('Constructed following ProcessInfoCollection from histogram list:')
  #print(PIC)

  # get valid processes and compare to arguments
  processes = args.processes
  if doallprocesses: processes = PIC.plist
  else:
    for p in processes:
      if p not in PIC.plist:
        raise Exception('ERROR: requested process {}'.format(p)
                        +' not found in the ProcessInfoCollection.')
  print('Extracted following valid process tags from input file:')
  for process in processes: print('  - '+process)
        
  # get valid systematics and compare to arguments
  shapesyslist = PIC.slist
  print('Extracted following relevant systematics from histogram file:')
  for systematic in shapesyslist: print('  - '+systematic)

  # loop over variables
  for var in varlist:
    print('Now running on variable {}...'.format(var.name))

    # get variable properties
    variablename = var.name
    variablemode = 'single'
    binlabels = None
    labelsize = None
    canvaswidth = None
    canvasheight = None
    p1legendbox = None
    p1legendncols = None
    labelangle = None
    if isinstance(var,DoubleHistogramVariable): variablemode = 'double'
    if variablemode=='single':
      xaxtitle = var.axtitle
      unit = var.unit
      if( var.iscategorical and var.xlabels is not None ):
        binlabels = var.xlabels
        labelsize = 15
    elif variablemode=='double':
      xaxtitle = var.primary.axtitle
      unit = var.primary.unit
      primarybinlabels = var.primary.getbinlabels()
      secondarybinlabels = var.secondary.getbinlabels(extended=True)
      binlabels = (primarybinlabels, secondarybinlabels)
      labelsize = 15
      labelangle = 45
      canvaswidth = 900
      p1legendbox = [0.45, 0.7, 0.95, 0.9]
      p1legendncols = 4

    # extra histogram selection for overlapping variable names
    othervarnames = [v.name for v in varlist if v.name!=variablename]
    thishistnames = lt.subselect_strings(histnames, 
                      mustcontainall=[variablename],
                      maynotcontainone=['_{}_'.format(el) for el in othervarnames])[1]

    # make a ProcessCollection for this variable
    splittag = args.region+'_'+variablename
    PIC = ProcessInfoCollection.fromhistlist( thishistnames, splittag, datatag=args.datatag )

    # manage systematics
    if( not args.rawsystematics and not args.dummysystematics ):
      _ = remove_systematics_default( PIC, year=args.year )
      _ = add_systematics_default( PIC, year=args.year )
    if args.dummysystematics:
      _ = remove_systematics_all( PIC )
      _ = add_systematics_dummy( PIC )

    # make a ProcessCollection
    PC = ProcessCollection( PIC, args.inputfile, doclip=True )

    # get the nominal simulated histograms
    simhists = []
    for process in PC.plist:
      simhists.append( PC.processes[process].hist )

    # now we have to rename the split histograms back to TTW0,...,TTW3
    if args.splitprocess is not None and args.splitvariable is not None:
      for hist in simhists:
        oldtitle = hist.GetTitle()
        lastchar = oldtitle[-1]
        if( lastchar.isdigit() ):
          hist.SetName(args.splitprocess + lastchar)
          hist.SetTitle(args.splitprocess + lastchar)

    # modify histogram titles
    for hist in simhists:
      title = hist.GetTitle()
      if title in processdict.keys():
        hist.SetTitle(processdict[title])

    # get the uncertainty histogram
    mcsysthist = PC.get_systematics_rss()

    # get data histogram
    datahistname = '{}_{}_{}_nominal'.format(args.datatag,args.region,variablename)
    if not datahistname in thishistnames:
      print('WARNING: no data histogram found.')
      datahist = PC.get_nominal()
      args.unblind = False
    else:
      f = ROOT.TFile.Open(args.inputfile,'read')
      datahist = f.Get(datahistname)
      datahist.SetDirectory(0) 
      f.Close()

    # blind data histogram
    if not args.unblind:
      for i in range(0,datahist.GetNbinsX()+2):
        datahist.SetBinContent(i, 0)
        datahist.SetBinError(i, 0)

    # set plot properties
    if( xaxtitle is not None and unit is not None ):
      xaxtitle += ' ({})'.format(unit)
    yaxtitle = 'Number of events'
    outfile = os.path.join(args.outputdir, variablename)
    lumimap = {'run2':137600, '2016':36300, '2017':41500, '2018':59700,
                    '2016PreVFP':19520, '2016PostVFP':16810 }
    if not args.year in lumimap.keys():
      print('WARNING: year {} not recognized,'.format(args.year)
            +' will not write lumi header.')
    lumi = lumimap.get(args.year,None)
    colormap = colors.getcolormap(style=args.colormap)
    if args.splitvariable is not None and args.splitprocess is not None and variablemode=='double':
        for key, value in colormap.items():
            if key[-1].isdigit():
                if int(key[-1])>0:
                    colormap[key+ args.splitvariable.replace('_','')] = value 
 
    extrainfos = []
    extrainfos.append( args.year )
    extrainfos.append( regionname )
    if args.splitvariable is not None and args.splitprocess is not None:
        extrainfos.append( args.splitprocess + " split on PL " + args.splitvariable )

    # for double histogram variables,
    # make a labelmap for better legends
    labelmap = None
    if( variablemode=='double' ):
      labelmap = {}
      # first 'regular' case where secondary variable is the split variable
      if args.splitvariable is None:
          splitvariable = var.secondary.name.strip('_')
          sbl_short = var.secondary.getbinlabels()
          for hist in simhists:
	      oldtitle = hist.GetTitle()
              newtitle = oldtitle[:]
              splitchar = ''
              # first case: names of the form TTW1
              if oldtitle[-1].isdigit():
                  splitchar = oldtitle[-1]
                  newtitle = newtitle[:-1]
              # second case: names of the form TTW1nMuons
              if oldtitle.endswith(splitvariable):
                  newtitle = newtitle.replace(splitvariable,'')
                  splitchar = newtitle[-1]
                  newtitle = newtitle[:-1]
                  # reset histogram title to correspond to first case (needed for correct colors)
                  oldtitle = newtitle+splitchar
                  hist.SetTitle(oldtitle)
	      if( splitchar.isdigit() ):
	          plbin = int(splitchar)
                  appendix = ''
	          if( plbin==0 ): appendix = '(o.a.)'
                  elif( plbin-1 < len(sbl_short) ): appendix = '({})'.format(sbl_short[plbin-1])
                  else: appendix = ''
                  if len(appendix)>0: newtitle = newtitle+' '+appendix
                  # also automatically add this process to the list of signals
                  signals.append(oldtitle)
              labelmap[oldtitle] = newtitle
      # now 'special' case where split variable can be different
      if( args.splitvariable is not None and args.splitprocess is not None ):  
          if args.splitvariable != var.secondary.name:
              print('WARNING: labels for case split variable != secondary variable not yet implemented.')
              continue
          signals = [x + args.splitvariable.replace('_','') for x in args.signals]
          signals.append(args.splitprocess+"0")
          signals.append(args.splitprocess+"1")
    
    # modify output file name as needed
    if( args.splitvariable is not None and args.splitprocess is not None ):
      outfile += '_split_' + args.splitvariable

    # make the plot
    hp.plotdatavsmc(outfile, datahist, simhists,
	    mcsysthist=mcsysthist, 
	    xaxtitle=xaxtitle,
	    yaxtitle=yaxtitle,
	    colormap=colormap,
            labelmap=labelmap,
            signals=args.signals,
            extrainfos=extrainfos,
	    lumi=lumi, extracmstext=args.extracmstext,
            binlabels=binlabels, labelsize=labelsize,
            labelangle=labelangle,
            canvaswidth=canvaswidth, canvasheight=canvasheight,
            p1legendbox=p1legendbox,
            p1legendncols=p1legendncols )

    if args.dolog:
      # make plot in log scale
      outfile = os.path.join(args.outputdir, variablename)+'_log'
      hp.plotdatavsmc(outfile, datahist, simhists,
            mcsysthist=mcsysthist,
            xaxtitle=xaxtitle,
            yaxtitle=yaxtitle,
            colormap=colormap,
            labelmap=labelmap,
            signals=args.signals,
            extrainfos=extrainfos,
            lumi=lumi, extracmstext=args.extracmstext,
            binlabels=binlabels, labelsize=labelsize,
            labelangle=labelangle,
            canvaswidth=canvaswidth, canvasheight=canvasheight,
            p1legendbox=p1legendbox,
            p1legendncols=p1legendncols,
            yaxlog=True )
