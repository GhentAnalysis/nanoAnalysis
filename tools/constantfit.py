#######################################################
# Translation of Tools/src/ConstantFit.cc into python #
#######################################################

import sys
import os
import ROOT
import math


def valueIsBad( val ):
  ### internal helper function to check failed fits
  if math.isnan(val): return True
  if math.isinf(val): return True  
  return False

def fitConstant( hist, xmin=None, xmax=None ):
  ### fit a constant to a TH1D

  # initialize result
  res = {'value': 0, 'uncertainty': 0, 'normalizedchi2': 0, 'success': False}

  # check if provided range is valid
  hnbins = hist.GetNbinsX()
  hxmin = hist.GetBinLowEdge(1)
  hxmax = hist.GetBinLowEdge(hnbins) + hist.GetBinWidth(hnbins)
  if xmin is None: xmin = hxmin
  if xmax is None: xmax = hxmax
  if( xmin < hxmin or xmax > hxmax ):
    print( "ERROR in fitConstant: given range falls outside of the given histogram's range." )
    return res

  # check for empty data
  if( hist.GetSumOfWeights() < 1e-6 ):
    print( "ERROR in fitConstant: histogram is empty." )
    return res

  # define the constant function to be fitted
  constfunc = ROOT.TF1( "constFunc", "[0]", xmin, xmax )

  # fit the function twice (hack to improve the result in root), in quiet mode ('Q')
  hist.Fit( "constFunc", "Q", "", xmin, xmax );
  hist.Fit( "constFunc", "Q", "", xmin, xmax );

  value = constfunc.GetParameter( 0 )
  uncertainty = constfunc.GetParError( 0 )
  normalizedchi2 = constfunc.GetChisquare() / constfunc.GetNDF()

  # check if fit was successful
  fitfailed = ( valueIsBad( value ) or valueIsBad( uncertainty ) or valueIsBad( normalizedchi2 ) )
  if fitfailed:
    print( "ERROR in fitConstant: fit failed." )
    return res

  # return the fit result
  res['value'] = value
  res['uncertainty'] = uncertainty
  res['normalizedchi2'] = normalizedchi2
  res['success'] = True
  return res
