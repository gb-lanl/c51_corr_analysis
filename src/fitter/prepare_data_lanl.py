import numpy as np
import gvar as gv
import re 
import pandas as pd 
import sys
import copy
# import tables as h5
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

from IPython import embed
sys.path.insert(0, '/home/gbradley/nucleon_elastic_FF')
from nucleon_elastic_ff.data.h5io import get_dsets 
from nucleon_elastic_ff.data.parsing import parse_t_info, parse_file_info 

from nucleon_elastic_ff.data.scripts.concat import concatenate
from nucleon_elastic_ff.data.scripts.concat import concat_dsets

''' 
Convert raw correlator data into agreeable format for fitter. 

Args: 


Returns:
    Coalesced_Dataset object: 
        corr_gv : dictionary of correlated data
'''

def coalesce_data(corr_raw, skip_prelim=False,fold=False,nt=None):

    corr_binned = raw_to_binned(
        corr_raw
    )

    corr_gv = Coalesced_Dataset(
        corr_binned,
        skip_prelim = skip_prelim,
        nt = nt
    )

    return corr_gv




def normalize_ff(curr,mom,m_snk):
    normalize = np.float(1)
    return normalize


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
    cnf_abbr = [files.split(".ama.h5")[0] for files in dirs]
    cnf_abbr = [cnf.replace('-','.') for cnf in cnf_abbr]
    cfg_abbr_sorted = np.sort(cnf_abbr)
    embed()
    
    # data_file_list = os.path.realpath(dirs)
    data_file_list = list()
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            data_file_list.append(os.path.abspath(os.path.join(dirpath, f)))
    sorted_files = np.sort(data_file_list)
    # embed()
    file = data_file_list[0]

    # cfg_pattern = "(E7-0_[0-9]+).ama.h5"
    # regex = re.compile(cfg_pattern)
    # df = pd.DataFrame([regex.match(f).groupdict() for f in data_file_list if regex.match(f)])
    # embed()



    string = '3pt_tsep21/NUCL_D_MIXED_NONREL_l0_g0/src10.0_snk10.0/qz+0_qy+0_qx+0/E7.a_1716/AMA'
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
    pattern_2pt = [
        "2pt/",
        "pion/",
        "src10.0_snk10.0/",
        "pion",]


    
    # for n in range(len(pattern_2pt)):
    #     pattern = "".join(pattern_2pt[:n+1])
    #     for i in range(len(data_file_list)):
    #         match = re.match(pattern_2pt, data_file_list[i])
    #         if not match:
    #             print(pattern_2pt)
    #             break

    # if match:
    #     concat_pattern = match.groupdict()
    # print(concat_pattern)
    # print(parse_t_info('3pt_tsep21/NUCL_D_MIXED_NONREL_l0_g0/src10.0_snk10.0/qz+0_qy+0_qx+0/E7.a_1716/AMA'))

    # concatenation_pattern = dict()
    # concatenation_pattern['2pt/pion/src10.0_snk10.0/pion/'] = 'pion'
    # concatenation_pattern['2pt/proton/src10.0_snk10.0/proton'] = 'proton'
    # concatenate(directory, concat_pattern)
    out_file = os.path.join(os.getcwd(),"conctat_out.h5")
    # print(out_file)
    # concat_dsets(data_file_list[0:3], out_file, overwrite=True)

#     string = (
#     "3pt_tsep12/NUCL_D_MIXED_NONREL_l0_g0/src5.0_snk5.0/qz+0_qy+0_qx+0/C13.b_5682/AMA"
# )

### Save data into txt files for unew analysis
    # newdir = os.path.join('unew_files',labels[ensemble]+'_pbp')
    # os.makedirs(newdir,exist_ok=True)
    # filename = os.path.join(newdir,'runs_')
    # for i,s in enumerate(streams[ensemble]):
    #     dset_mdt[s].to_csv(filename+str(i)+'.dat',sep=' ',header=False,index=False,columns=['pbp_l','pbp_s','pbp_c'],mode='w')

    for f_i in range(0,3):
        with h5py.File(sorted_files[f_i], 'r') as h5f:
            dsets = get_dsets(h5f,load_dsets=True)
            dsets.to_csv()
            # embed()
            # print(dsets.keys())
            
            # concat_dsets(data_file_list, out_file,overwrite=True)


        # embed()
    # pion = {}
    for cfg in cfg_abbr_sorted[0:2]:
        
        pion   = dsets['2pt/pion/src10.0_snk10.0/pion/'+cfg+'/AMA']
        # proton = dsets['2pt/proton/src10.0_snk10.0/proton/E7.a_1716/AMA']
    print(pion)
        # embed()
    def get_real(data):
        pion_real = []
        for i in range(len(data)):
            pion_real.append(data[i][0])
        return pion_real

    raw = {}
    raw['pion'] = get_real(pion)
    print(raw['pion'])
    print(gv.dataset.avg_data(raw['pion']))
    mean = gv.mean(gv.dataset.avg_data(raw['pion']))
    # cov = correct_covariance(data, **kwargs)
    print(gv.gvar(mean))
        # data = gv.gvar(pion_real)
    # print(np.mean(pion_real),np.std(pion_real))

    def bs_corr(corr,Nbs,Mbs,seed=None):
        corr_bs = np.zeros(tuple([Nbs]) + corr.shape[1:],dtype=corr.dtype)
        np.random.seed(seed) # if None - it does not seed - I checked 14 May 2013
        # make bs_lst of shape (Nbs,Mbs)
        bs_lst = np.random.randint(0,corr.shape[0],(Nbs,Mbs))
        # use bs_lst to make corr_bs entries
        for bs in range(Nbs):
            corr_bs[bs] = corr[bs_lst[bs]].mean(axis=0)
        return corr_bs

    # bs_corr(pion_real, Nbs=20, Mbs=20,seed='a071m170')

    
    # to_gvar(data)

    def effective_mass(data):
        """
        Computes the effective mass analytically using the following formula
        
        meff = ArcCosh( (C(t+1)+ C(t-1)) / C(t) )
        This method correctly accounts for contributions both from forward- and
        backward-propagating states. 
        """
        cosh_m = (data[2:] + data[:-2]) / (2.0 * data[1:-1])
        meff = np.zeros(len(cosh_m), dtype=gv._gvarcore.GVar)
        # The domain of arccosh is [1, Infinity).
        # Set entries outside of the domain to nans.
        domain = (cosh_m > 1)
        meff[domain] = np.arccosh(cosh_m[domain])
        meff[~domain] = gv.gvar(np.nan)
        return meff

    effective_mass(data)
            # corrs_gv = {}
            # corrs_gv = gv.dataset.avg_data(pion_real)
            # print(corrs_gv)
            

        
       
        # threept = dsets['3pt_tsep21/NUCL_D_MIXED_NONREL_l0_g0/src10.0_snk10.0/qz+0_qy+0_qx+0/E7.a_1716/AMA']
        # corrs = gv.BufferDict()
        # corrs['pion'] = dsets['2pt/pion/src10.0_snk10.0/pion/E7.a_1716/AMA']
        # corrs['proton'] = dsets['2pt/proton/src10.0_snk10.0/proton/E7.a_1716/AMA']
        # corrs_gv = {}
        # for k in corrs:
        #     corrs_gv[k] = gv.dataset.avg_data(corrs[k])
        # #     else:
        #         for corr in corr_dict:
        #             corrs_corr = {k: v for k, v in corrs.items() if corr in k}
        #             tmp_gv = gv.dataset.avg_data(corrs_corr)
        #             for k in tmp_gv:
        #                 corrs_gv[k] = tmp_gv[k]
        # else:
        #     corrs_gv = gv.dataset.avg_data(corrs)


        # print(dsets)
        # for key, dset in dsets.items():
        #     match = re.search(pattern, string)
        #     if match:
        #         info = match.groupdict()

if __name__== '__main__':
    main()