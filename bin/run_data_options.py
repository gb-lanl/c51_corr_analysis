import sys 
import argparse
import argcomplete

from fitter.prepare_data_lanl import * 
from fitter.load_data import * 
#from utils import *

class Format(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def main():
    # formatting for help messages 
    formatter = lambda form: argparse.RawTextHelpFormatter(form) 

    parser = argparse.ArgumentParser(description= "Compute averages of provided data"
    "with several statistical options to choose from")
    parser.add_argument('-filename(s)', nargs='*',default=[sys.stdin])
    
    

    parser.add_argument('-bs', '--bootstrap')
    parser.add_argument('-jn', '--jackknife')
    parser.add_argument('-sdev', '--standard_deviation')

    parser.add_argument('-columns', '--data_columns')
    parser.add_argument('-jn_nbs', '--jackknife_blocks')

    args = parser.parse_args()
    # add path to the input file and load it
    # if no sys arg given for filename but input file instead:
    # sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
    # fp = importlib.import_module(
    #     args.fit_params.split('/')[-1].split('.py')[0])


    """ read in pure """ 
    """ naive read in """
    if args.bootstrap:
    
    if args.jackknife:

    if args.standard_deviation:

if __name__ == '__main__':
    main()