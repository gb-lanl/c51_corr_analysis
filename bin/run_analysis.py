import sys
import lsqfit
import gvar as gv
import importlib
import h5py
import copy
import os
import pathlib
import random
import matplotlib.pyplot as plt
import argparse
import numpy as np
import argparse
import pandas as pd
np.seterr(invalid='ignore') 

# lanl_analysis_suite libs.modules
from fitter import prepare_data_lanl as ld
from fitter import corr_functions as cf
from fitter import plotting as visualize 
import utilities



def main():
        print(*sys.argv)
        #     reparse_argv()
        less_indent_formatter = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=40)
        parser = argparse.ArgumentParser(description = 
        " This is the master executable for analysis of \n"
        " lqcd correlator and form factor data. The fallback default arguments are contained in /tests/input_file/defaults.py. \n"
        "The logical procedure is as follows:" 
        "Features have been added with modularity/user flexibility in mind. However, with this type of analysis,"
        " The human-provided input remains necessary. ")
        
        parser.add_argument('fit_params',    help='input file to specify fit')
        parser.add_argument('-b', '--block', default=1, type=int,
                        help=            'specify bin/blocking size in terms of saved configs')
        parser.add_argument('--states',      nargs='+',
                        help=            'specify states to fit?')
        parser.add_argument('--plot-file', default = None, dest = 'res_filename',
                help = "Plot the correlator from a file instead of performing a fit. \n"
                " You can pass a fitmass... file here")
        parser.add_argument('--plot-start', action = 'store_true', dest = 'plot_start',
                help = "Do not perform a fit. Instead generate a plot with the start parameters.\n"
                " Has to be passed along with --start-params")


        # parser.add_argument('filename', help = "The filename containing the data")
        # parser.add_argument('n_states override', help = "override given n_states in input file here \n")
        # parser.add_argument('--Nt', 
        #         help = "dont compute Nt from data, use this one. WARNING: correlator \n"
        # "will not be symmetrized")
        # parser.add_argument('--log-level', default = "INFO",
        #         help = "log level options: WARN, INFO, PROGRESS, DETAILS, DEBUG, NONE \n")
        # parser.add_argument('n_states override', help = "override given n_states in input file here \n")

        # parser.add_argument('channel', help='Select channels to fit')
        # parser.add_argument('run_src', help= "Run fit for C_2pt at the src?")
        # parser.add_argument('run_snk', help= "Run fit for C_2pt at the snk?")
        # parser.add_argument('run_ratio', help= "Run fit for C_3pt/C_2pt?")
        # parser.add_argument('run_direct', help= "Run fit for spectral decomposition (not using C_2pt)?")
        # parser.add_argument('summary', help="print summary of fit results in current shell?") #TODO to txt, save data pickle





        args = parser.parse_args()
        # TODO set up buffers
        # Buffer = 
        # buffers = []
        # Channel = namedtuple("Channel", ['abbr','name'])
        # channels = []

        ''' parse provided input/label file '''
        sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
        fp = importlib.import_module(
                args.fit_params.split('/')[-1].split('.py')[0])
        
        # block data
        bl = args.block
        if 'block' in dir(fp) and fp.block != 1:
                bl = fp.block
        if args.block != 1: # allow cl override
                bl = args.block

        if args.states:
                states = args.states
        else:
                states = fp.fit_states

        # corr_raw = 

        # ''' process correlator data '''
        # corr_gv = fitter.prepare_data_lanl.coalesce_data(corr_raw)

        # ''' data handling options ''' 



        h5f = h5py.File('src/fitter/proton_all.h5','r')
        h5f_3pt = h5py.File('src/fitter/3pt.h5','r')
        # h5f_sp = h5py.File('proton_SP.h5','r')
        string = '2pt/proton/src10.0_snk10.0/proton/AMA'
        string_sp = '2pt/proton_SP/src10.0_snk10.0/proton/AMA'
        string_tsep13_d = '3pt_tsep13/NUCL_D_l0_g4/AMA'
        string_tsep13_u = '3pt_tsep13/NUCL_U_l0_g4/AMA'
        c3pt_data = {'A3': {}, 'V4' : {} }
        c2pt_data = {}
        # sinks = [""]
        # for g_type in ['g7','g11','g13','g14']:
        #         for T in [13,15,17,19,21]:
        #                 string_tsep13_d = '3pt_tsep13/NUCL_D_l0_g4/AMA'
        #                 string_tsep13_u = '3pt_tsep13/NUCL_U_l0_g4/AMA'
        #                 string_tsep15 = '3pt_tsep15/NUCL_D_l0_g4/AMA'
        #                 string_tsep17 = '3pt_tsep17/NUCL_D_l0_g4/AMA'
        #                 string_tsep19 = '3pt_tsep19/NUCL_D_l0_g4/AMA'
        #                 string_tsep21 = '3pt_tsep21/NUCL_D_l0_g4/AMA'

        to_array = lambda f, path : pd.DataFrame(f[path][:]).values #real for v4, imag else 
        temp_data = (to_array(h5f_3pt,string_tsep13_u) - to_array(h5f_3pt,string_tsep13_d))
        c3pt_data[int(13)] = temp_data
        print(len(c3pt_data[int(13)]))
        ifil = h5f[string][:]
        ifil_sp = h5f[string_sp][:]
        # ifil_13 = h5f_3pt[string_tsep13][:]
        
        # only need real part for 2pt correlators
        _ifil = pd.DataFrame(ifil)['re'] 
        _ifil_sp = pd.DataFrame(ifil_sp)['re'] 
        corrs = {}
        corrs['proton'] = _ifil.to_numpy()
        corrs['proton_SP'] = _ifil_sp.to_numpy()
        # corrs[int(13)] = _ifil_13.to_numpy()
        # corrs.columns
        Ncfg = corrs['proton'].shape[0]
        bl = 18
        Ncfg_3pt = 8232 #c3pt_data[int(13)].shape[0]
        bl_ = 294
        
        ydata = {}
        ydata_3pt = {}
        ydata['proton_SS'] = np.reshape(corrs['proton'], (Ncfg//bl, bl))
        ydata['proton_SP'] = np.reshape(corrs['proton_SP'], (Ncfg//bl, bl))
        ydata_3pt[int(13)] = np.reshape(c3pt_data[int(13)], (Ncfg_3pt//bl_, bl_))
        print(ydata_3pt[int(13)])
        # ydata[int(15)] = 
        # ydata[int(17)] = 
        # ydata[int(19)] =
        # ydata[int(21)] = 

        # print(ydata.keys())
        ydata_out = ld.bs_to_gvar(data=ydata, corr='proton_SS',bs_N=100) #576?

        x, priors = ld.prepare_xyp(states, fp, ydata_out)
        # print(x,priors)
        c2pt = {}
        c2_src = cf.C_2pt(tag='proton_SS', ydata=ydata_out['proton_SS'],p=priors)
        c2_snk = cf.C_2pt(tag='proton_SP', ydata=ydata_out['proton_SP'],p=priors)

        print(c2_src.meff())

        c2pt['SS'] = c2_src
        c2pt['PS'] = c2_snk
        print(c2_src.times)
        c2_avg = {tag: c2pt[tag].avg() for tag in c2pt}
        print(c2_avg)
        tfit = {}
        for tag in c2pt:
                # print(tag)
                tfit[tag] =  np.arange(c2pt[tag].times.tmin,c2pt[tag].times.tmax +1)
        print(tfit)

        # print(c2_src.plot_corr())
        print(visualize.get_nucleon_effective_mass(ydata_out))
        visualize.plot_effective_mass(ydata_out)
        
        c3 = cf.C_3pt('A3', ydata_3pt=ydata_3pt[int(13)],t_ins=range(16),T=[13,15,17])

        """ time containers for correlators """
        nts = [c2_src.times.nt,c2_snk.times.nt, c3.times.nt]
        print(nts)

        tmaxes = [c2_src.times.tmax, c2_snk.times.tmax, c3.times.tmax]
        tdata = np.arange(0, min(tmaxes))

        """ Estimate ground-state mass associated with source operator."""
        m_src = c2_src.mass

        """ Estimate ground-state mass associated with sink operator."""
        m_snk = c2_snk.mass
        if m_src is not None:
                c2_src.set_mass(m_src)
        if m_snk is not None:
                c2_snk.set_mass(m_snk)

   
        # print(c3.avg(m_src=c2_src.meff(), m_snk=c2_snk.meff()))
        # print(ydata_out)
        # for k,v in ydata_out.items():
        #         print(k,v)
        # print(states)
        
        


        # lg.set_log_level(args.log_level)
        # args = vars(args)
        # del args['log_level']

        # '''Run fits sequentially'''

        # if args.run_src:
        #         analysis.run_src(n=args.nstates)

        # if args.run_snk:
        #         analysis.run_snk(n=args.nstates)

        
        # if args.run_ratio:
        #         analysis.run_ratio(
        #                 n=args.nstates-1, 
        #                 tmin_src = c2_src.times.tmin,  
        #                 tmin_snk = c2_snk.times.tmin,
        #                 t_iter = c3.times.nt)
        #         # if args.axial:
        #         #         analysis
                
        #         # if args.vector:
                
        #         # if args.scalar:
                
        #         # if args.tensor:

        # if args.run_direct:
        #         analysis.run_direct(nstates=args.nstates)
        
        # if args.summary:
        #         analysis.print_summary()





        











if __name__ == '__main__':
    main()