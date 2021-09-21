import argparse
import os 
import sys


# Running together with the script:
# ================================== #

#     #!/bin/bash

#     #SBATCH --job-name=ggp
#     #SBATCH --cpus-per-task=1
#     #SBATCH --mem-per-cpu=4G

#     #SBATCH --time=23:00:00
#     #SBATCH --qos=1day

#     #SBATCH --output=std.out
#     #SBATCH --mail-type=END,FAIL,TIME_LIMIT

#     ${COMMAND}


# Example
# ---------------------------------- #
# python3 run_all_Theo.py -d ../data_theo/ -b ../data_theo/parameters_GM.txt ../data_theo/parameters_SPM.txt -c ../data_theo/csv_config.txt -m -t 1e-20

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


def run_command(arg, dryrun, iscluster):
    # run on cluster
    if iscluster:
        com = "sbatch  --export=COMMAND='" + arg +"'" + " submit_ggp_run.sl"
        if dryrun:  
            print(com)
        else:
            print(arg)
            os.system(com)

    # run locally
    else:
        if dryrun:  
            print(arg)
        else:
            os.system(arg)        

def get_arg_list(args):
    s = ''
    for a in args:
        s += ' ' + a
    return s + ' '


def get_parameter_files(file):
    out_dir = file[:-4] + '_out/'
    entries = os.listdir(out_dir)
    final_files = []
    for e in entries:
        if e.endswith('parameter_file.txt'):
            final_files.append(os.path.join(out_dir, e))
    
    # sort accoring to segment number (files only differ at that number)
    print(len(final_files))
    if len(final_files)>1:
        final_files = sorted(final_files)

    return final_files

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

    parser.add_argument('-b',
                        dest='parameters' ,
                        help='parameter file',
                        nargs='+',
                        required=True)

    parser.add_argument('-c',
                        dest='csv_config' ,
                        help='csv_config file',
                        nargs='+',
                        required=True)

    parser.add_argument('-space',
                        dest='space' ,
                        help='space',
                        default='log')

    parser.add_argument('-t',
                        dest="tol",
                        help="tolerance",
                        required=False,
                        default="1e-7")
    
    parser.add_argument('--dryrun', action='store_true')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('-m', action='store_true')
    parser.add_argument('-p', action='store_true')

    args = parser.parse_args()

    # ======================================== #
    # ======================================== #

    if args.local:
        config =    {"bin": "./../bin/gfp_gaussian",
                     "iscluster": False}
    else:
        config =    {"bin": "./../bin/gfp_gaussian",
                     "iscluster": True}


    input_files = get_input_files(args.dir)

    for infile in input_files:
        ggp_arg =   " -c " + get_arg_list(args.csv_config) + \
                    " -t " + args.tol + \
                    " -space " + args.space + \
                    " -i " + infile 
        if args.m:
            ggp_arg += ' -m '
        if args.p:
            ggp_arg += ' -p '

        # deal with parameter files
        if args.parameters[0] == "infer":
            parameter_files = get_parameter_files(infile) 
        else:
            parameter_files = args.parameters

        ggp_arg +=  " -b " + get_arg_list(parameter_files)
        
        # ============ run! ============ #
        run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"])


# ================================================================================ #
if __name__ == "__main__":
    main()