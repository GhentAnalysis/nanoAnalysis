###################################################
# Concrete implementation of Run-II UL reweighter #
###################################################

import sys
import os
import numpy as np
import awkward as ak
from pathlib import Path
# import all reweighters
sys.path.append(str(Path(__file__).parents[1]))
from combinedreweighter import CombinedReweighter
from muonrecoreweighter import MuonRecoReweighter
from muonidreweighter import MuonIDReweighter
from electronidreweighter import ElectronIDReweighter
from electronrecoreweighter import ElectronRecoReweighter
from pileupreweighter import PileupReweighter
from prefirereweighter import PrefireReweighter
from customreweighter import CustomReweighter
from btagreweighter import BTagReweighter
from psreweighter import PSReweighter
from pdfreweighter import PDFReweighter
from scalereweighter import ScaleReweighter


def njets_evaluator(events, **kwargs):
    # initialize
    nominal = np.ones(len(events))
    up = np.ones(len(events))
    down = np.ones(len(events))
    # get jets
    if 'jet_mask' not in kwargs.keys():
        msg = 'ERROR: reweighter.weight was called without'
        msg += ' required argument "jet_mask".'
        raise Exception(msg)
    jets = events.Jet[kwargs['jet_mask']]
    njets = ak.num(jets)
    # 10% uncertainty per jet up to 50%
    for i in [1,2,3,4,5]:
      up = np.where(njets>=i, 1+i/10., up)
      down = np.where(njets>=i, 1-i/10., down)
    return {'nominal': nominal, 'up': up, 'down': down}

def nbjets_evaluator(events, **kwargs):
    # initialize
    nominal = np.ones(len(events))
    up = np.ones(len(events))
    down = np.ones(len(events))
    # get jets
    if 'bjet_mask' not in kwargs.keys():
        msg = 'ERROR: reweighter.weight was called without'
        msg += ' required argument "bjet_mask".'
        raise Exception(msg)
    bjets = events.Jet[kwargs['bjet_mask']]
    nbjets = ak.num(bjets)
    # 10% uncertainty for 0 or 1 b-jets, 40% for 2+ b-jets
    up = np.where(nbjets>=0, 1.1, up)
    down = np.where(nbjets>=0, 0.9, down)
    up = np.where(nbjets>=2, 1.4, up)
    down = np.where(nbjets>=2, 0.6, down)
    return {'nominal': nominal, 'up': up, 'down': down}

def get_run2ul_reweighter( 
      year,
      sampleweights,
      dobtagnormalize=True ):
    ### get a correctly configured combined reweighter
    # input arguments:
    # - year: data taking year (in str format)
    # - sampleweights: an object of type SampleWeights for the current sample
    # - dobtagnormalize: set the b-tagging reweighter to use normalization

    # initializations
    reweighter = CombinedReweighter()
    weightdir = os.path.join(os.path.dirname(__file__), '../data')

    # electron reco reweighter
    wfile = os.path.join(weightdir, 'electronreco', 
              'electronreco_sf_{}.json'.format(year))
    reweighter.add_reweighter('electronreco', ElectronRecoReweighter(wfile, year))

    # muon reco reweighter
    wfile = os.path.join(weightdir, 'muonreco', 
              'muonreco_sf_{}.json'.format(year))
    #reweighter.add_reweighter('muonreco', MuonRecoReweighter(wfile))

    # electron id reweighter
    wfile = os.path.join(weightdir, 'leptonid',
              'leptonMVAUL_SF_electrons_Tight_{}.root'.format(year))
    reweighter.add_reweighter('electronid', ElectronIDReweighter(wfile))

    # muon id reweighter
    wfile = os.path.join(weightdir, 'leptonid',
              'leptonMVAUL_SF_muons_Medium_{}.root'.format(year))
    reweighter.add_reweighter('muonid', MuonIDReweighter(wfile, 'Medium'))

    # b-tag reweighter
    wfile = os.path.join(weightdir, 'btagging', 'btagging_{}.json'.format(year))
    reweighter.add_reweighter('btagging', BTagReweighter(wfile, normalize=dobtagnormalize))

    # pileup reweighter
    wfile = os.path.join(weightdir, 'pileup', 'puWeights_{}.json'.format(year))
    reweighter.add_reweighter('pileup', PileupReweighter(wfile, year))

    # prefire reweighter
    reweighter.add_reweighter('prefire', PrefireReweighter(prefiretype='all'))

    # njets and nbjets reweighter
    reweighter.add_reweighter('njets', CustomReweighter(njets_evaluator))
    reweighter.add_reweighter('nbjets', CustomReweighter(nbjets_evaluator))

    # parton shower reweighter
    reweighter.add_reweighter('ps', PSReweighter())

    # scale reweighters
    scalereweighter = ScaleReweighter(sampleweights, rtype='acceptance')
    reweighter.add_reweighter('scaleacceptance', scalereweighter)
    scalereweighter = ScaleReweighter(sampleweights, rtype='norm')
    reweighter.add_reweighter('scalenorm', scalereweighter)

    # pdf reweighters
    pdfreweighter = PDFReweighter(sampleweights, rtype='acceptance')
    reweighter.add_reweighter('pdfacceptance', pdfreweighter)
    pdfreweighter = PDFReweighter(sampleweights, rtype='norm')
    reweighter.add_reweighter('pdfnorm', pdfreweighter)

    return reweighter
