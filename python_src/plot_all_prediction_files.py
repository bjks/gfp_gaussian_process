import argparse
import os 
import sys

from read_ggp_run import *


def mk_missing_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory) 
    return directory

def get_input_files(directory, keyword=None):
    entries = os.listdir(directory)
    final_files = []
    if keyword == None:
        for e in entries:
            if e.endswith(".csv"):
                final_files.append(os.path.join(directory,e))
    else:
        for e in entries:
            if e.endswith(".csv") and keyword in e:
                final_files.append(os.path.join(directory,e))   

    return final_files


def get_prediction_file(directory):
    entries = os.listdir(directory)
    final_files = []
    for e in entries:
        if e.endswith('prediction.csv'):
            final_files.append( os.path.join(directory, e) )
        # if e.endswith('backward.csv'):
        #     final_files.append( os.path.join(directory, e) )
        # if e.endswith('forward.csv'):
        #     final_files.append( os.path.join(directory, e) )
    return final_files

########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Plot predictions')

    parser.add_argument('-d',
                        dest='dir',
                        help='dir',
                        required=True)

    parser.add_argument('-skip',
                        dest='skip',
                        help='skip rows in input file',
                        default=25,
                        type=int,
                        required=False)
    
    parser.add_argument('-r',
                        dest='range',
                        help='range (start, stop, step)',
                        nargs='+',
                        type=int,
                        default=[0, 100, 1],
                        required=False)
    
    parser.add_argument('-o',
                        dest='out',
                        help='Out dir',
                        required=False)


    args = parser.parse_args()
    # ======================================== #
    # ======================================== #

    input_files = get_prediction_file(args.dir)

    mk_missing_dir(args.out)

    for infile in input_files:
        print(infile, 'cells:', args.range[0], args.range[1], args.range[2])

        out_file = os.path.join(args.out, infile.split("/")[-1][:-4]) + ".png"
        plot_predictions(infile, start=args.range[0], stop=args.range[1], step=args.range[2], 
                    time_unit=("h", 1), skip_row=args.skip, xlim=[None, None], outfile=out_file)


# ================================================================================ #
if __name__ == "__main__":
    main()