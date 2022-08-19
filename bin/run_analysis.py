import sys 
from fitter import *


def main():
    print(*sys.argv)
    reparse_argv()
    less_indent_formatter = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=40)
    parser = argparse.ArgumentParser(description = " This is the master executable for analysis of \n"
    " lqcd correlator fits. The fallback default arguments are contained in  The logical procedure is as follows: \n"
    " 1. \n" 
    "Features have been added with modularity/user flexibility in mind. However, with this type of analysis,"
    " The human- provided input remains necessary. ")
    

    parser.add_argument('--plot-file', default = None, dest = 'res_filename',
            help = "Plot the correlator from a file instead of performing a fit. \n"
            " You can pass a fitmass... file here")
    parser.add_argument('--plot-start', action = 'store_true', dest = 'plot_start',
            help = "Do not perform a fit. Instead generate a plot with the start parameters.\n"
            " Has to be passed along with --start-params")


    parser.add_argument('filename', help = "The filename containing the data")
    pasrser.add_argument('n_states override', help = "override given n_states in input file here \n")
    pasrser.add_argument('--Nt', 
            help = "dont compute Nt from data, use this one. WARNING: correlator \n"
    "will not be symmetrized")
    pasrser.add_argument('--log-level', default = "INFO",
            help = "log level options: WARN, INFO, PROGRESS, DETAILS, DEBUG, NONE \n")
    pasrser.add_argument('n_states override', help = "override given n_states in input file here \n")
    





    args = parser.parse_args()

    lg.set_log_level(args.log_level)
    args = vars(args)
    del args['log_level']



if __name__ == '__main__':
    main()