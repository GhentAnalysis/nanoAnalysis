#################################
# Tools for usage with argparse #
#################################

import sys
import os

def path_or_none( p ):
    ### return None if p is None, else take absolute path
    if p is None: return None
    else: return os.path.abspath(p)

def path_or_das( p ):
    #### return raw value if p points to a file on DAS, else take absolute path
    if p.startswith('root://'): return p
    else: return os.path.abspath(p)
