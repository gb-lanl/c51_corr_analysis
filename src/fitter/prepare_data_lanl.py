import numpy as np
import gvar as gv
import re 
import pandas as pd 
import sys
import copy
import tables as h5
import h5py
import os 
import time
import re
import argparse
import importlib
import corrfitter
import collections
import numpy as np
import lsqfit

sys.path.insert(0, '/home/gbradley/nucleon_elastic_FF')
from nucleon_elastic_ff.data.h5io import get_dsets 

def main():
    parser = argparse.ArgumentParser(
        description='Perform analysis of two-point correlation function')
    parser.add_argument('fit_params',    help='input file to specify fit')
    parser.add_argument('--fit',         default=True, action='store_true',
                        help=            'do fit? [%(default)s]')
    parser.add_argument('--prior_override',default=True, action='store_true',
                        help=            'override generated priors with priors from input file? [%(default)s]')
    parser.add_argument('-b', '--block', default=1, type=int,
                        help=            'specify bin/blocking size in terms of saved configs')
    parser.add_argument('--uncorr_corrs', default=False, action='store_true',
                        help=            'uncorrelate different correlation functions? [%(default)s]')
    parser.add_argument('--uncorr_all',  default=False, action='store_true',
                        help=            'uncorrelate all snk,src for each correlation function? [%(default)s]')
    parser.add_argument('--states',      nargs='+',
                        help=            'specify states to fit?')    
    parser.add_argument('--svdcut',      type=float, help='add svdcut to fit')
    parser.add_argument('--svd_test',    default=True, action='store_false',
                        help=            'perform gvar svd_diagnosis? [%(default)s]')
    parser.add_argument('--svd_nbs',     type=int, default=50, help='number of BS samples for estimating SVD cut [%(default)s]')                    

    args = parser.parse_args()
    # if args.save_figs and not os.path.exists('figures'):
    #     os.makedirs('figures')
    # print(args)
    # add path to the input file and load it
    sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
    fp = importlib.import_module(
        args.fit_params.split('/')[-1].split('.py')[0])

    directory = fp.data_directory
    N_cnf = len([name for name in os.listdir(directory) if os.path.isfile(name)])

    dirs = os.listdir( directory )
    cnf_abbr = [files.split(".ama.h5",0) for files in dirs]

    # data_file_list = os.path.realpath(dirs)
    data_file_list = list()
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            data_file_list.append(os.path.abspath(os.path.join(dirpath, f)))
    file = data_file_list[0]

    string = (
    "3pt_tsep12/NUCL_D_MIXED_NONREL_l0_g0/src5.0_snk5.0/qz+0_qy+0_qx+0/C13.b_5682/AMA"
)

    patterns = [
        "3pt",
        "_tsep(?P<tsep>[0-9]|[0-9]+)",  # must match `_tsep` and stores the following numbers (any length)
        "/NUCL_(?P<quark>U|D)",  # Store U or D in quark
        "_MIXED_NONREL",  # Not sure if this changes. Not stored for now
        "_l(?P<l>[0-9]+)",  # action parameters?
        "_g(?P<g>[0-15]+)",
        "/src(?P<src>[0-9\.]+)",  # Stores numbers + . to store decimals. Must escape .
        "_snk(?P<snk>[0-9\.]+)",  # Stores numbers + . to store decimals. Must escape .
        "/qz(?P<qz>[\+\-0-9]+)", 
        "_qy(?P<qy>[\+\-0-9]+)", 
        "_qx(?P<qx>[\+\-0-9]+)", 
        
    ]
    with h5py.File(file, "r") as h5f:
        dsets = get_dsets(h5f)
        # print(dsets)
        for key, dset in dsets.items():
            match = re.search(pattern, string)
            if match:
                info = match.groupdict()

if __name__== '__main__':
    main()