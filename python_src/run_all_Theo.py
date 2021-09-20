import argparse
import os 
import sys


# slurm_script = \
# """
# #!/bin/bash

# #SBATCH --job-name=ggp
# #SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=4G

# #SBATCH --time=23:00:00
# #SBATCH --qos=1day

# #SBATCH --output=std.out
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT


# """

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

def get_new_parameter_files(input_file):
    base = input_file.split('/')[-1][:-4]
    p0 = input_file[:-4] + '_out/' + base + '_segment0_f012345678910_b_parameter_file.txt'
    p1 = input_file[:-4] + '_out/' + base + '_segment1_f012345678910_b_parameter_file.txt'
    return p0, p1



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


    for file in get_input_files(args.dir):
        ggp_arg =   " -b " + get_arg_list(args.parameters) + \
                    " -c " + get_arg_list(args.csv_config) + \
                    " -t " + args.tol + \
                    " -space " + args.space + \
                    " -i " + file 
        if args.m:
            ggp_arg += ' -m '
        if args.p:
            ggp_arg += ' -p '

        command = config["bin"] + ' ' + ggp_arg

        run_command(command, args.dryrun, config["iscluster"])


# ================================================================================ #
if __name__ == "__main__":
    main()