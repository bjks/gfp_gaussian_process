import argparse
import os 
import sys

from read_ggp_run import *
import glob

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


def get_prediction_files(path):
    if os.path.isfile(path):
        return [path]
    candidates = sorted(glob.glob(path + '/**/*prediction.csv', recursive=True))
    # ignore directories that end with '.czi' and only select files!
    return [c for c in candidates if os.path.isfile(c)]

########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Plot predictions')

    parser.add_argument('-d',
                        dest='dir',
                        help='Directory(-ies) that will be searched for prediction files (recursively)',
                        nargs='+',
                        type=str,
                        required=True)

    parser.add_argument('-time_unit',
                        dest='time_unit',
                        help='Time unit and rescaling of the time eg t/60, default [min, 1]',
                        nargs='+',
                        default=["min", "1"],
                        type=str,
                        required=False)
    
    parser.add_argument('-r',
                        dest='range',
                        help='Range for cells that will be plotted (start, stop, step)',
                        nargs='+',
                        type=int,
                        default=[None, None, 1],
                        required=False)
    
    parser.add_argument('-o',
                        dest='out',
                        help='Output directory, default: same as prediction file',
                        nargs='+',
                        type=str,
                        default=[],
                        required=False)

    parser.add_argument('--replot', help="Replot plots for files, that already exist", action='store_true')


    args = parser.parse_args()
    # ======================================== #
    # ======================================== #
    for i, directory in enumerate(args.dir):
        input_files = get_prediction_files(directory)

        if len(args.out)>i:
            mk_missing_dir(args.out[i])

        for infile in input_files:
            print(infile, 'cells:', args.range[0], args.range[1], args.range[2])

            if len(args.out)>i:
                out_file = os.path.join(args.out[i], infile.split("/")[-1][:-4]) + ".png"
            else:
                out_file = infile[:-4] + ".png"

            if args.replot or not os.path.exists(out_file):
                plot_predictions(infile, 
                                start=args.range[0], stop=args.range[1], step=args.range[2], 
                                time_unit=(args.time_unit[0], float(args.time_unit[1])), 
                                skip_row = header_lines(infile, until="cell_id"), 
                                xlim=[None, None], 
                                outfile=out_file)


# ================================================================================ #
if __name__ == "__main__":
    main()