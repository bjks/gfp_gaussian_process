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


def get_prediction_file(file):
    out_dir = file[:-4] + '_out/'
    entries = os.listdir(out_dir)
    final_files = []
    for e in entries:
        if e.endswith('prediction.csv'):
            return os.path.join(out_dir, e)
    return None

def get_plot_file(file):
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
        description='Run gfp_gaussian for all files in dir')

    parser.add_argument('-d',
                        dest='dir',
                        help='dir',
                        required=True)
    
    parser.add_argument('--dryrun', action='store_true')
    parser.add_argument('--local', action='store_true')


    args = parser.parse_args()

    # ======================================== #
    # ======================================== #


    input_files = get_input_files(args.dir)

    for infile in input_files:
        print(infile)
        prediction_file = get_prediction_file(infile)

        if prediction_file!=None:
            plot_predictions(prediction_file, start=10, stop=100, step=1, 
                    time_unit=("min", 60), skip_row=25, xlim=[None, None], outfile=get_plot_file(infile))


# ================================================================================ #
if __name__ == "__main__":
    main()