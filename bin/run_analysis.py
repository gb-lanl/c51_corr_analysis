import sys 
from src import fitter
from tests import *
import argparse




def main():
        print(*sys.argv)
        #     reparse_argv()
        less_indent_formatter = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=40)
        parser = argparse.ArgumentParser(description = " This is the master executable for analysis of \n"
        " lqcd correlator and form factor data. The fallback default arguments are contained in /tests/input_file/defaults.py. \n"
        "The logical procedure is as follows:" 
        "Features have been added with modularity/user flexibility in mind. However, with this type of analysis,"
        " The human-provided input remains necessary. ")
        

        parser.add_argument('--plot-file', default = None, dest = 'res_filename',
                help = "Plot the correlator from a file instead of performing a fit. \n"
                " You can pass a fitmass... file here")
        parser.add_argument('--plot-start', action = 'store_true', dest = 'plot_start',
                help = "Do not perform a fit. Instead generate a plot with the start parameters.\n"
                " Has to be passed along with --start-params")


        parser.add_argument('filename', help = "The filename containing the data")
        parser.add_argument('n_states override', help = "override given n_states in input file here \n")
        parser.add_argument('--Nt', 
                help = "dont compute Nt from data, use this one. WARNING: correlator \n"
        "will not be symmetrized")
        parser.add_argument('--log-level', default = "INFO",
                help = "log level options: WARN, INFO, PROGRESS, DETAILS, DEBUG, NONE \n")
        parser.add_argument('n_states override', help = "override given n_states in input file here \n")

        parser.add_argument('channel', help='Select channels to fit')
        parser.add_argument('run_src', help= "Run fit for C_2pt at the src?")
        parser.add_argument('run_snk', help= "Run fit for C_2pt at the snk?")
        parser.add_argument('run_ratio', help= "Run fit for C_3pt/C_2pt?")
        parser.add_argument('run_direct', help= "Run fit for spectral decomposition (not using C_2pt)?")
        parser.add_argument('summary', help="print summary of fit results in current shell?") #TODO to txt, save data pickle





        args = parser.parse_args()
        ''' parse provided input/label file '''
        # corr_raw = 

        ''' process correlator data '''
        corr_gv = fitter.prepare_data_lanl.coalesce_data(corr_raw)

        ''' data handling options ''' 



        lg.set_log_level(args.log_level)
        args = vars(args)
        del args['log_level']

        '''Run fits sequentially'''

        if args.run_src:
                analysis.run_src(n=args.nstates)

        if args.run_snk:
                analysis.run_snk(n=args.nstates)

        ''' nts = [c2_src.times.nt,
           c2_snk.times.nt,
           c3.times.nt]

        tmaxes = [c2_src.times.tmax,
                c2_snk.times.tmax,
                c3.times.tmax]
        tdata = np.arange(0, min(tmaxes)) '''
        if args.run_ratio:
                analysis.run_ratio(
                        n=args.nstates-1, 
                        tmin_src = c2_src.times.tmin,  
                        tmin_snk = c2_snk.times.tmin,
                        t_iter = c3.times.nt)
                # if args.axial:
                #         analysis
                
                # if args.vector:
                
                # if args.scalar:
                
                # if args.tensor:

        if args.run_direct:
                analysis.run_direct(nstates=args.nstates)
        
        if args.summary:
                analysis.print_summary()





        











if __name__ == '__main__':
    main()