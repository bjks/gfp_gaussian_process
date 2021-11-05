import argparse
import os 
import sys

from read_ggp_run import *

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
    return final_files

def get_plot_file(files):
    out_dir = file[:-4] + '_out/'
    entries = os.listdir(out_dir)
    final_files = []
    for e in entries:
        if e.endswith('prediction.csv'):
            return os.path.join(out_dir, e[:-4]+'.png')
    return None

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


    args = parser.parse_args()
    # ======================================== #
    # ======================================== #


    input_files = get_input_files(args.dir)

    for infile in input_files:
        print(infile, 'cells:', args.range[0], args.range[1], args.range[2])
        prediction_file = get_prediction_file(infile)

        if prediction_file!=None:
            plot_predictions(prediction_file, start=args.range[0], stop=args.range[1], step=args.range[2], 
                    time_unit=("min", 60), skip_row=args.skip, xlim=[None, None], outfile=get_plot_file(infile))


# ================================================================================ #
if __name__ == "__main__":
    main()