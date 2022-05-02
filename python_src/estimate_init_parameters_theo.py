import argparse
import os 
import sys
import pandas as pd
import numpy as np

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


def estimate_lambda(cells):
    lambdas = np.array([])
    for cell in cells:
        dx = np.diff(np.log(cell.length))
        dt = np.diff(cell.time)
        
        lambdas = np.append(lambdas, dx/dt)
    return np.max([np.mean(lambdas), 1e-6]), np.var(lambdas)

def estimate_q(cells, beta=0):
    qs = np.array([])
    for cell in cells:
        dg = np.diff(cell.gfp)
        dt = np.diff(cell.time)

        length_interp = cell.length[:-1] + np.diff(cell.length)/2
        gfp_interp = cell.gfp[:-1] + np.diff(cell.gfp)/2

        qs = np.append(qs, (dg/dt + beta*gfp_interp)/length_interp)

    if np.mean(qs) < 0:
        print('mean q', np.mean(qs), '!!!!!')

    return np.max([np.mean(qs), 0.01]), np.var(qs)

def estimate_var_x(cells, rel_err=0.05):
    x = []
    for cell in cells:
        x = np.append(x, np.log(cell.length))
    return (rel_err * np.mean(x))**2

def estimate_var_g(cells, rel_err=0.05):
    g = []
    for cell in cells:
        g = np.append(g, cell.gfp)
    return (rel_err * np.mean(g))**2

def estimate_var_dx(cells, rel_dev=0.5):
    x = []
    for cell in cells:
        x = np.append(x, np.log(cell.length[0]))
    return (rel_dev * np.mean(x))**2

def estimate_var_dg(cells, rel_dev=0.5):
    g = []
    for cell in cells:
        g = np.append(g, cell.gfp[0])
    return (rel_dev * np.mean(g))**2

def write_params2file(params, filname):
    with open(filname, 'w') as fout:
        fout.write("# Automatically estimated parameter for initialzing MLE search\n")
        for param in params:
            if params[param][0] == "bound":
                fout.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format(param, params[param][1],
                                                                 params[param][1]/2.,
                                                                 params[param][2],
                                                                 params[param][3]) )
            elif params[param][0] == "free":
                fout.write("{:s} = {:.2E}, {:.2E}\n".format(param, params[param][1], 
                                                                    params[param][1]/2.))
            elif params[param][0] == "fixed":
                fout.write("{:s} = {:.2E}\n".format(param, params[param][1]))


def get_paramter_file(directory, infile, seg_idx):
    entries = os.listdir(directory)
    for e in entries:
        if infile.split("/")[-1][:-4] in e  and "_parameter_file.txt" in e and "segment" + str(seg_idx):
            return os.path.join(directory, e)

def get_var_g(filename):
    with open(filename, 'r') as fout:
        for line in fout:
            if line.startswith("var_g"):
                return line.split("=")[-1]



def look_for_output_file(out_dir, infile, outype):
    entries = os.listdir(out_dir)
    for e in entries:
        core = infile.split('/')[-1][:-4]
        if core in e and outype+".csv" in e:
                return os.path.join(out_dir, e)
    return None
########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Estimate init parameters from raw input files.\n Example: python3 estimate_init_parameters.py -d ../../experimental_data/data_theo_20211213/experimental_data_GFP -s GM.1_parameters_GFP SPM.1_parameters_GFP')

    parser.add_argument('-d',
                        dest='dir',
                        help='Directory with input files',
                        required=True)

    
    parser.add_argument('-rescale_time',
                        dest='rescale_time',
                        help='rescale_time',
                        type=float,
                        default=3600,
                        required=False)

    parser.add_argument('-suffix',
                        dest='segments',
                        help='suffix for segments',
                        nargs='+',
                        type=str,
                        default=["parameters"],
                        required=False)

    parser.add_argument('-look_up_prediction',
                        dest='look_up_prediction',
                        help='Directory where the predition file might be, if not existing write parameter files',
                        default=None,
                        required=False)

    parser.add_argument('-gamma_lambda',
                        dest='gamma_lambda',
                        help='guessed gamma_lambda',
                        type=float,
                        default=[5, 0.01, 100],
                        required=False)    

    parser.add_argument('-gamma_q',
                        dest='gamma_q',
                        help='guessed gamma_q',
                        nargs='+',
                        type=float,
                        default=[5, 0.01, 100],
                        required=False)                 

    parser.add_argument('--free_beta', help="Set beta to 'free'", action='store_true')

    args = parser.parse_args()

    # ======================================== #
    # ======================================== #

    input_files = get_input_files(args.dir)
    print(len(input_files), "input files in", args.dir)

    for infile in input_files:
        if args.look_up_prediction!=None:
            p_file = look_for_output_file(args.look_up_prediction, infile, "prediction")
            if p_file!=None:
                print("Prediction file exists, do not rewrite parameter files")
                continue

        if args.free_beta:
            type_beta = "free"
        else:
            type_beta = "fixed"

        if "GFP" in infile:
            beta_bounds = [[type_beta, 0.01], # GM
                            [type_beta, 0.005]] # SP
        elif "RFP" in infile:
            beta_bounds = [[type_beta, 0.13],
                            [type_beta, 0.018]]  
        elif "YFP" in infile:
            beta_bounds = [[type_beta, 0.3],
                            [type_beta, 0.08]]  

        data = pd.read_csv(infile, skiprows=0)
        data.loc[: , "time"] = data["time_sec"].to_numpy() / args.rescale_time

        cells_data_segs = []
        params = {}

        for seg, prefix in enumerate(args.segments):
            print(infile, seg)
            if 'segment' in data.columns:
                segment_data = data[data["segment"] == seg]
            else:
                segment_data = data

            cells_data = df2raw_cells(segment_data, 
                                    time="time", 
                                    length="length_um", 
                                    gfp="fluo_ampl_ch_1", 
                                    cell_id="cell_ID", 
                                    parent_id="parent_ID")

            # cells_data = df2raw_cells(segment_data, 
            #                         time="time", 
            #                         length="length_um", 
            #                         gfp="gfp_nb", 
            #                         cell_id="sub_cell", 
            #                         parent_id="sub_parent")

            # OUs
            mean_lambda, variance_lambda = estimate_lambda(cells_data)

            params["mean_lambda"] = ["free", mean_lambda]
            params["gamma_lambda"] = ["bound"] + args.gamma_lambda
            params["var_lambda"] = ["free", variance_lambda/(2*args.gamma_lambda[0])]


            mean_q, variance_q = estimate_q(cells_data, beta_bounds[seg][1]) 

            params["mean_q"] = ["free", mean_q]
            params["gamma_q"] = ["bound"] + args.gamma_q
            params["var_q"] = ["free", variance_q/(2*args.gamma_q[0])]
            
            params["beta"] =  beta_bounds[seg]

            # measurment noise
            params["var_x"] = ["free", estimate_var_x(cells_data, rel_err=0.01)]
            params["var_g"] = ["free", 1]

            # cell division
            params["var_dx"] = ["free", estimate_var_dx(cells_data, rel_dev=0.5)]
            params["var_dg"] = ["free", estimate_var_dg(cells_data, rel_dev=0.5)]

            # write to file named as expected by "run_all_Theo.py"
            parameter_file = infile[:-4] + '_'+ prefix + ".txt"
            write_params2file(params, parameter_file)

if __name__ == "__main__":
    main()