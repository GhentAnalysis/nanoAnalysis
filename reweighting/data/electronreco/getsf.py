#################################################
# Small utility script to get the scale factors #
#################################################
# see https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON

import os
import sys

if __name__=='__main__':

    years = {
      '2016preVFP': '2016PreVFP',
      '2016postVFP': '2016PostVFP',
      '2017': '2017',
      '2018': '2018'
    }
    topdir = '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM'

    for key, val in years.items():
        fname = os.path.join(topdir, key+'_UL', 'electron.json.gz')
        if not os.path.exists(fname):
            raise Exception('ERROR: file {} does not exist.'.format(fname))
        newfname = 'electronreco_sf_{}.json.gz'.format(val)
        cmd = 'cp {} {}'.format(fname, newfname)
        os.system(cmd)
        cmd = 'gzip -d {}'.format(newfname)
        os.system(cmd)
