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
                        help='Directory with input files',
                        required=True)

    parser.add_argument('-o',
                        dest='out',
                        help='Do not use the default out dir but use this directory instead',
                        default=None,
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
        try:
            outfile = get_plot_file(infile, args.out)

            # plot_raw_data_input_file(infile, 
            #                 time="time_sec", 
            #                 length="length_um", 
            #                 gfp="fluo_ampl_ch_1", 
            #                 cell_id="cell_ID", 
            #                 parent_id="parent_ID", 
            #                 start=args.range[0], stop=args.range[1], step=args.range[2], 
            #                 time_unit=("min", 60), outfile=outfile, 
            #                 scatter=False)
            # print(infile, '->', outfile, 'cells:', args.range[0], args.range[1], args.range[2])

            plot_raw_data_input_file(infile, 
                            time="time_sec", 
                            length="length_um", 
                            gfp="gfp_nb", 
                            cell_id="sub_cell", 
                            parent_id="sub_parent", 
                            start=args.range[0], stop=args.range[1], step=args.range[2], 
                            time_unit=("min", 60), outfile=outfile, 
                            scatter=False)
            print(infile, '->', outfile, 'cells:', args.range[0], args.range[1], args.range[2])




        except:
            print( "\n ----------> ", infile, "FAILED! <----------\n" )
# ================================================================================ #
if __name__ == "__main__":
    main()