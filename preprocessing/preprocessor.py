################################################################
# Convenience class grouping commonly used preprocessing steps #
################################################################

# imports
import sys
import os
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))
import preprocessing.leptonvariables as lepvars
import preprocessing.leptongenvariables as lepgenvars
import preprocessing.topleptonmva as topmva
import preprocessing.triggervariables as triggervars


class PreProcessor(object):
    ### convenience class grouping commonly used preprocessing steps
    
    def __init__( self ):
        ### initializer
        # empty for now, maybe extend later
        pass

    def process( self, events,
                 leptonvariables=None,
                 leptongenvariables=None,
                 topmvavariable=None, topmvaversion=None,
                 dotriggers=False ):
        ### do preprocessing of a set of events
        # - use leptonvariables = ['all'] to calculate all lepton variables.
        # - use leptongenvariables = ['all'] to calculate all lepton generator variables.
        # - calculating the TOP lepton MVA score is disabled by default;
        #   use topmvavariable = 'mvaTOP' (or another name) 
        #   and topmvaversion = 'ULv1' (or another version)
        #   to enable it.

        # do checks
        if( topmvavariable is not None and 'year' not in events.metadata.keys() ):
            msg = 'ERROR in PreProcessor.process:'
            msg += ' to calculate TOP MVA scores, the year must be set in the events metadata.'
            raise Exception(msg)

        # add additional lepton variables
        if( leptonvariables is not None and len(leptonvariables)>0 ):
            lepvars.add_electron_variables(events, variables=leptonvariables)
            lepvars.add_muon_variables(events, variables=leptonvariables)

        # add additional gen-info variables
        if( leptongenvariables is not None and len(leptongenvariables)>0 ):
            # extra check on data type
            if 'GenPart' not in events.fields:
                msg = 'ERROR in PreProcessor.process:'
                msg += ' generator variables are requested, but sample appears to be data.'
                raise Exception(msg)
            lepgenvars.add_electron_gen_variables(events, variables=leptongenvariables)
            lepgenvars.add_muon_gen_variables(events, variables=leptongenvariables)

        # add TOP lepton mva scores
        if topmvavariable is not None:
            reader = topmva.TopLeptonMvaReader(events.metadata['year'], topmvaversion, verbose=True)
            reader.set_scores(events, name=topmvavariable)

        # add aggregated triggers
        if dotriggers:
            triggervars.add_trigger_variables(events, year=events.metadata['year'])
