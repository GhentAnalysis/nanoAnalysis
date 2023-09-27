##################################################
# Classes for samples and collections of samples #
##################################################

import sys
import os


def sample_is_2016(sample):
    if( '16PreVFP' in sample or '16PostVFP' in sample ): return False
    return (
      'Run2016' in sample # data
    )

def sample_is_2016PreVFP(sample):
    return (
      '16PreVFP' in sample # general
      or 'RunIISummer20UL16NanoAODAPV' in sample # UL simulation
      or 'Run2016B' in sample or 'Run2016C' in sample # data
      or 'Run2016D' in sample or 'Run2016E' in sample # data
      or ('Run2016F' in sample and 'HIPM' in sample) # data
    )

def sample_is_2016PostVFP(sample):
    if sample_is_2016PreVFP(sample): return False
    return (
      '16PostVFP' in sample # general
      or 'RunIISummer20UL16' in sample # UL simulation
      or 'Run2016F' in sample or 'Run2017G' in sample # data
      or 'Run2016H' in sample # data
    )

def sample_is_2017(sample):
    return (
      'RunIISummer20UL17' in sample # UL simulation 
      or 'Run2017' in sample # data
    )

def sample_is_2018(sample):
    return (
      'RunIISummer20UL18' in sample # UL simulation 
      or 'Run2018' in sample # data
    )

def year_from_sample_name(sample):
    ### internal helper function to extract year from a sample file name
    # input arguments:
    # - sample: sample file name in str format
    # returns:
    # year in str format
    years = []
    if sample_is_2016PreVFP(sample): years.append('2016PreVFP')
    if sample_is_2016PostVFP(sample): years.append('2016PostVFP')
    if sample_is_2017(sample): years.append('2017')
    if sample_is_2018(sample): years.append('2018')
    if len(years)==0:
        msg = 'ERROR in samples.year_from_sample_name:'
        msg += ' could not find any valid year for sample {}.'.format(sample)
        raise Exception(msg)
    if len(years)>1:
        msg = 'ERROR in samples.year_from_sample_name:'
        msg += ' found multiple valid years for sample {}: {}.'.format(sample, years)
        raise Exception(msg)
    return years[0]

def dtype_from_sample_name(sample):
    ### internal helper function to extract sim or data from a sample file name
    # input arguments:
    # - sample: sample file name in str format
    # returns:
    # data type in str format ('sim' or 'data')
    dtypes = []
    if( 'RunII' in sample ): dtypes.append('sim')
    elif( 'Run2016' in sample
        or 'Run2017' in sample
        or 'Run2018' in sample ): dtypes.append('data')
    if len(dtypes)==0:
        msg = 'ERROR in samples.dtype_from_sample_name:'
        msg += ' could not find any valid dtype for sample {}.'.format(sample)
        raise Exception(msg)
    if len(dtypes)>1:
        msg = 'ERROR in samples.dtype_from_sample_name:'
        msg += ' found multiple valid dtypes for sample {}: {}.'.format(sample, dtypes)
        raise Exception(msg)
    return dtypes[0]


class Sample(object):

    def __init__( self ):
        self.name = None
        self.process = None
        self.path = None
        self.xsec = 0.
        self.year = None
        self.dtype = None

    def read_from_line( self, line, sampledir=None,
        doyear=True, dodtype=True, **kwargs ):
        ### read sample properties from a sample list line.
        # input arguments:
        # - line: string representing a line from a sample list.
        #         expected format: <process> <file> <xsec>
        #         notes: - cross-section is optional.
        #                - file can be either the full absolute path to the file,
        #                  or a file (base) name or relative path w.r.t. sampledir.
        #                - in extension, file can also be the name of a dataset on DAS.
        # - sampledir: if specified, set the full sample path and check if it exists.
        #              should be set to None if the samples are absolute local paths,
        #              or if they represent dataset names on DAS.
        # - doyear: if True, determine sample year from file name.
        # - dodtype: if True, determine data type (sim or data) from file name.
        # - kwargs: passed down to Sample.set_path

        # split the line by spaces
        line = line.strip(' ').split(' ')
        # remove extra spaces
        line = [el for el in line if el!='']
        # first extract the tag (short name) of the process
        self.process = line[0]
        # now extract sample file
        self.name = line[1].rstrip('\n')
        # finally extract cross-section
        if len(line)>2:
            xsstr = line[2].rstrip('\n')
            try: self.xsec = float(xsstr)
            except: print('WARNING in Sample.read_from_line:'
                            +' found incompatible cross-section "{}";'.format(xsstr)
                            +' using zero as default.')
        # set the path for this sample
        if sampledir is not None: self.set_path( sampledir, **kwargs )
        # get the year for this sample
        if doyear: self.year = year_from_sample_name(self.name)
        # get the data type for this sample
        if dodtype: self.dtype = dtype_from_sample_name(self.name)

    def set_path( self, sampledir, suppress_exception=False ):
        ### set the path attribute
        path = os.path.join(sampledir, self.name)
        if not os.path.exists(path):
            if suppress_exception: return
            raise Exception('ERROR in Sample.set_path:'
                            +' path {} does not exist.'.format(path))
        self.path = path

    def __str__( self ):
        res = 'Sample( name: {}'.format(self.name)
        res += ', process: {}'.format(self.process)
        res += ', xsec: {}'.format(self.xsec)
        res += ', path: {}'.format(self.path)
        res += ', year: {}'.format(self.year)
        res += ', dtype: {})'.format(self.dtype)
        return res
        

class SampleCollection(object):

    def __init__( self ):
        self.samples = []

    def read_from_file( self, samplelistpath, **kwargs ):
        self.read_from_files( [samplelistpath], **kwargs )

    def read_from_files( self, samplelistpaths, **kwargs ):
        ### read sample collection from a sample list file.
        # input arguments:
        # - samplelistpaths: list of sample list paths
        # - kwargs: passed down to Sample.read_from_line
        for samplelist in samplelistpaths:
            with open(samplelist) as f:
                for line in f:
                    # strip
                    line = line.strip(' \t\n')
                    # ignore blank or commented lines
                    if(len(line)<=1): continue
                    if(line[0] == '#'): continue
                    # make a sample
                    sample = Sample()
                    sample.read_from_line(line, **kwargs)
                    # add the sample to the list
                    self.samples.append(sample)

    def get_samples( self ):
        ### return the sample collection
        return self.samples

    def number( self ):
        ### return number of samples
        return len(self.get_samples())

    def check_paths( self, return_list=False ):
        ### check if all samples have an existing path assigned
        # input arguments:
        # - return list: if False, return a bool (True if all good, False otherwise);
        #                if True, return a list of samples with no valid path.
        missing_samples = []
        for sample in self.samples:
            if( sample.path is None or not os.path.exists(sample.path) ):
                missing_samples.append(sample)
        if return_list: return missing_samples
        else: return (len(missing_samples)==0)

    def __str__( self ):
        return '\n'.join(['{}'.format(s) for s in self.samples])
