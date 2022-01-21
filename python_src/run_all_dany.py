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


def get_parameter_file(file, suffix = "_parameters"):
    parameter_file = file[:-4] + suffix + '.txt'
    return parameter_file


########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Run gfp_gaussian for all files in dir (written for Dany's data sets)")

    parser.add_argument('-d',
                        dest='dir',
                        help='directory with input files',
                        required=True)
    
    parser.add_argument('-o',
                        dest='out' ,
                        help='Output dir',
                        default=None)

    parser.add_argument('-c',
                        dest='csv_config' ,
                        help='Csv config file',
                        required=True)

    parser.add_argument('-space',
                        dest='space' ,
                        help='Space, log(default) or linear',
                        default='log')

    parser.add_argument('-t',
                        dest="tol",
                        help="Tolerance of maximization",
                        default="1e-7")
    
    parser.add_argument('--dryrun', help="Shows what will be done", action='store_true')
    parser.add_argument('--local', help="Do not submit job, but run directly", action='store_true')

    parser.add_argument('-m', help="Run maximization", action='store_true')
    parser.add_argument('-p', help="Run prediction", action='store_true')
    parser.add_argument('-j', help="Run prediction", action='store_true')

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
        ggp_arg =   " -c "      + args.csv_config + \
                    " -t "      + args.tol + \
                    " -space "  + args.space + \
                    " -i "      + infile 
        if args.m:
            ggp_arg += ' -m '
        if args.p:
            ggp_arg += ' -p '
        if args.j:
            ggp_arg += ' -j '

        parameter_file = get_parameter_file(infile) 

        ggp_arg +=  " -b " + parameter_file

        if args.out != None:
            ggp_arg +=  " -o " + args.out

        # ============ run! ============ #
        run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"])


# ================================================================================ #
if __name__ == "__main__":
    main()