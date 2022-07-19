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

    return sorted(final_files)



def get_plot_file(file, out):
    if out==None:
        out_dir = file[:-4] + '_out/'
        return out_dir + "raw_traces.png"
    else:
        base = file[:-4].split('/')[-1]
        if not os.path.exists(out):
            os.mkdir(out)
        return os.path.join(out, base) + "raw_traces.png"

########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Plot raw data')

    parser.add_argument('-d',
                        dest='dir',
                        help='Directory that will be searched for input files',
                        required=True)

    parser.add_argument('-o',
                        dest='out',
                        help='Output directory, default: same as prediction file',
                        default=None,
                        required=False)

    parser.add_argument('-r',
                        dest='range',
                        help='Range for cells that will be plotted (start, stop, step)',
                        nargs='+',
                        type=int,
                        default=[None, None, 1],
                        required=False)

    parser.add_argument('-time_unit',
                        dest='time_unit',
                        help='Time unit and rescaling of the time eg t/60, default [min, 1]',
                        default=["min", "1"],
                        type=str,
                        required=False)

    ### Reading the csv                    
    parser.add_argument('-time_col',
                        dest='time_col',
                        help='Column, default: time_col',
                        default="time_sec",
                        required=False)
    
    parser.add_argument('-length_col',
                        dest='length_col',
                        help='Column, default: length_um',
                        default="length_um",
                        required=False)

    parser.add_argument('-gfp_col',
                        dest='gfp_col',
                        help='Column, default: gfp_nb',
                        default="gfp_nb",
                        required=False)

    parser.add_argument('-cell_id',
                        dest='cell_id',
                        help='Column, default: cell_ID',
                        default="cell_ID",
                        required=False)

    parser.add_argument('-parent_id',
                        dest='parent_id',
                        help='Column, default: parent_ID',
                        default="parent_ID",
                        required=False)



    args = parser.parse_args()

    # ======================================== #
    # ======================================== #

    input_files = get_input_files(args.dir)

    for infile in input_files:
        try:
            outfile = get_plot_file(infile, args.out)



            plot_raw_data_input_file(infile, 
                            time = args.time_col, 
                            length = args.length_col, 
                            gfp = args.gfp_col, 
                            cell_id = args.cell_id, 
                            parent_id = args.parent_id, 
                            start = args.range[0], stop=args.range[1], step=args.range[2], 
                            time_unit=(args.time_unit[0], float(args.time_unit[1])), 
                            outfile=outfile, 
                            scatter=False)

            print(infile, ' plot is saved in', outfile, ', cells:', args.range[0], args.range[1], args.range[2])

        except Exception as e:
            print("ERROR :", str(e), ";", infile, "FAILED\n")
# ================================================================================ #
if __name__ == "__main__":
    main()