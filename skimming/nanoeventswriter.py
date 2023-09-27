####################################################
# Class for writing an events array to a root file #
####################################################
# This class performs the opposite operation as 
# coffea.nanoevents.NanoEventsFactory.from_root (with the NanoAODSchema).
# Where the latter builds an events array from a NanoAOD root file,
# this class writes the array back to a compatible root file,
# potentially after event selection, dropping unused branches and/or adding custom branches.
# The resulting file can be read again with the NanoEventsFactory
# for further processing steps.

import sys
import os
import awkward as ak
import uproot
import ROOT


class NanoEventsWriter(object):

    def __init__( self ):
        ### initializer
        # empty for now, maybe extend later
        pass

    def write( self, events, outputfile, compressionlevel=9, drop=None ):
        ### write events to a root file
        # input arguments:
        # - events: a NanoEventsArray object
        # - outputfile: name of outputfile to write (must be .root)
        # - drop: tuple of the form (fields to drop, subfields to drop),
        #   where fieds to drop is a list (or None) of fields to discard,
        #   and subfields to drop is a dict (or None) of field names to subfields to discard.
        #   drop can also be a str identifier for a predefined set of fields and subfields to drop
        #   (see implementation below).
        # documentation: this method is based on:
        # https://uproot.readthedocs.io/en/latest/basic.html#writing-ttrees-to-a-file
        
        # get lists of branches to drop
        fieldstodrop = []
        subfieldstodrop = {}
        if drop is not None:
            if isinstance(drop, str):
                fieldstodrop, subfieldstodrop = self.get_fields_to_drop(drop)
            else:
                if drop[0] is not None: fieldstodrop = drop[0]
                if drop[1] is not None: subfieldstodrop = drop[1]

        # transform the events object to a tree dict
        tree = {}
        for field in sorted(events.fields):
            if field in fieldstodrop: continue
            if( len(events[field].fields)>0 ):
                subtree = {}
                for subfield in sorted(events[field].fields):
                    # skip global index mappings that are created by NanoEventsFactory,
                    # but that are not present in original nanoAOD file
                    if subfield.endswith('IdxG'): continue
                    # do a dummy assertion to make sure that the branch is materialized
                    # (if turned off, can give nasty and hard-to-debug errors upon writing).
                    # as an alternative to assert, one could also print the counts.
                    assert ak.count(events[field][subfield])>=0
                    subtree[subfield] = events[field][subfield]
                tree[field] = ak.zip(subtree)
            else:
                tree[field] = events[field]

        # write to file
        with uproot.recreate(outputfile, compression=uproot.LZMA(compressionlevel)) as f:
            f['Events'] = tree
        

    def get_fields_to_drop( self, dropid ):
        ### internal helper function for getting sets of fields to drop
        if dropid=='default':
            fieldstodrop = (['CaloMET','ChsMET','CorrT1METJet',
              'DeepMETResolutionTune','DeepMETResponseTune',
              'FatJet','FsrPhoton',
              'GenDressedLepton','GenIsolatedPhoton','GenJetAK8',
              'GenMET','GenVisTau','GenVtx',
              'HTXS','IsoTrack','L1','L1Reco','L1simulation',
              'LHEPart','LowPtElectron','OtherPV',
              'PuppiMET','RawMET','RawPuppiMET','SV',
              'SubGenJetAK8','SubJet','TkMET','boostedTau'])
            subfieldstodrop = None
            return (fieldstodrop, subfieldstodrop)
        # throw error if id was not recognized
        msg = 'ERROR in NanoEventsWriter.get_fields_to_drop:'
        msg += ' identifer {} was not recognized.'.format(dropid)
        raise Exception(msg)


class AuxiliaryTreeWriter(object):

    # note: easiest using simple PyROOT?
    # note: strange errors on a test file for the ParameterSets tree with unknown cause...
    #       ignore the MetaData and ParameterSets tree for now.

    def __init__( self ):
        ### initializer
        # empty for now, maybe extend later
        pass

    def write( self, inputfile, outputfile ):
        ### write auxiliary nanoAOD trees to a root file
        # input arguments:
        # - inputfile: name of inputfile to read (must be .root)
        # - outputfile: name of outputfile to write (must be .root)
        #   note: outputfile can be an existing root file,
        #         in which case the trees will be added to the file,
        #         without deleting its earlier contents.
        f = ROOT.TFile.Open(inputfile)
        lumitree = f.Get( "LuminosityBlocks" ).CloneTree()
        runtree = f.Get( "Runs" ).CloneTree()
        lumitree.SetDirectory(0)
        runtree.SetDirectory(0)
        f.Close()
        f = ROOT.TFile.Open(outputfile, 'update')
        lumitree.Write()
        runtree.Write()
        f.Close()
