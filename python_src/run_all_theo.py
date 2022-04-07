import argparse
import os 
import sys
from tkinter import E

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
# python3 run_all_Theo.py -d ../data_theo/ -b ../data_theo/parameters_GM.txt ../data_theo/parameters_SPM.txt --fallback -c ../data_theo/csv_config.txt -m -t 1e-20
# python3 run_all_Theo.py -d ../data_theo/ -b infer -c ../data_theo/csv_config_NOFILTER.txt -p --local

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


def run_command(arg, dryrun, iscluster, verbose):
    # run on cluster
    if iscluster:
        com = "sbatch  --export=COMMAND='" + arg +"'" + " submit_ggp_run.sl"
        if dryrun:
            if verbose:  
                print(arg)
        else:
            # print(arg)
            os.system(com)
    # run locally
    else:
        if dryrun:
            if verbose:  
                print(arg)
        else:
            os.system(arg)        

def get_arg_list(args):
    s = ''
    for a in args:
        s += ' ' + a
    return s + ' '


def get_parameter_files(file, out_dir = None):
    if out_dir == None:
        out_dir = file[:-4] + '_out/'
    entries = os.listdir(out_dir)
    parameter_files = []
    for e in entries:
        if file.split("/")[-1][:-4] in e and e.endswith('parameter_file.txt'):
            parameter_files.append(os.path.join(out_dir, e))
    
    # sort accoring to segment number (files only differ at that number)
    if len(parameter_files)>1:
        parameter_files = sorted(parameter_files)

    return parameter_files


def free_up_parameters(paramter_file, free_params):
    file_data = []
    with open(paramter_file,'r') as fin:
        for _, line in enumerate(fin):
            file_data.append(line)
    with open(paramter_file,'w') as fout:
        for line in file_data:
            if line.split("=")[0].strip() in free_params:
                param = line.split("=")[0].strip() 
                init = str(float(line.split("=")[-1].split(',')[0]))
                step = str(float( init )/2.)
                newline = param + ' = ' + init + ' , ' + step + "\n"
                fout.write(newline)
                if param.startswith("gamma_q"):
                    newline = param + ' = ' + init + ' , ' + step + ' , ' + '0.01' + ' , ' + '100' "\n"
            else:
                fout.write(line)
    print("Rewrote", paramter_file)



# def create_new_parameter_file_from_both(parameter_files):
#     segment0 = []
#     segment1 = []
#     if len(parameter_files) == 2:
#         with open(parameter_files[0],'r') as fin:
#             for _, line in enumerate(fin):
#                 if not line.startswith('#'):
#                     segment0.append(line)
#         with open(parameter_files[1],'r') as fin:
#             for _, line in enumerate(fin):
#                 if not line.startswith('#'):
#                     segment1.append(line)

#         output_file = parameter_files[0].replace("segment0", "segment2" )
#         with open(output_file, 'w') as fout:
#             fout.write("#Automatically generated file identical to segments1 apart from beta which from segments0\n")
#             for i,_ in enumerate(segment0):
#                 if segment0[i].startswith("beta"):
#                     fout.write(segment0[i])
#                 else:
#                     fout.write(segment1[i])
#         return parameter_files + [output_file]
#     else:
#         return parameter_files 

# def create_new_parameter_file(parameter_files):
#     segment1 = []
#     if len(parameter_files) == 2:
#         with open(parameter_files[1],'r') as fin:
#             for _, line in enumerate(fin):
#                 if not line.startswith('#'):
#                     segment1.append(line)

#         output_file = parameter_files[0].replace("segment0", "segment2" )
#         with open(output_file, 'w') as fout:
#             fout.write("#Automatically generated file identical to segments1 apart from beta\n")
#             for i,_ in enumerate(segment1):
#                 if segment1[i].startswith("beta"):
#                     line_splitted = segment1[i].split("=")
#                     fout.write(line_splitted[0] + ' = ' + str(float(line_splitted[1])*4) + '\n')
#                 else:
#                     fout.write(segment1[i])
#         return parameter_files + [output_file]
#     else:
#         return parameter_files 


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
        description="Run gfp_gaussian for all files in dir (written for Theo's data sets)")

    parser.add_argument('-d',
                        dest='dir',
                        help='directory with inpput files',
                        required=True)

    parser.add_argument('-b',
                        dest='parameters' ,
                        help='Parameter file(s)',
                        nargs='+',
                        default=[],
                        required=False)

    parser.add_argument('-free_up',
                        dest='free_up' ,
                        help='Free up parameters in SP',
                        nargs='+',
                        default=[],
                        required=False)


    parser.add_argument('-suffix',
                        dest='suffix' ,
                        help='Parameter suffixes',
                        nargs='+',
                        default=[],
                        required=False)

    parser.add_argument('-o',
                        dest='out' ,
                        help='Output dir (None)',
                        default=None)

    parser.add_argument('-c',
                        dest='csv_config' ,
                        help='Csv config file',
                        required=True)

    parser.add_argument('-space',
                        dest='space' ,
                        help='Space, log or linear (log)',
                        default='log')

    parser.add_argument('-t',
                        dest="tol",
                        help="Tolerance of maximization (1e-12)",
                        required=False,
                        default="1e-12")
    
    parser.add_argument('-noise',
                        dest="noise",
                        help="Noise model (scaled)",
                        required=False,
                        default="scaled")
    
    parser.add_argument('--dryrun', help="Shows what will be done", action='store_true')
    parser.add_argument('--verbose', help="Shows what will be done", action='store_true')

    parser.add_argument('--local', help="Do not submit job, but run directly", action='store_true')
    # parser.add_argument('--fallback', help="Indicate that the parameter files are fallbacks if there are no specific files", action='store_true')
    # parser.add_argument('--newparamfile', help="Creates third parameter file ('segment2') that is identical to the segment1 but contains the beta from segment0", action='store_true')
    parser.add_argument('--rerun', help="Rerun for files, that already have a prediction file", action='store_true')

    parser.add_argument('-m', help="Run maximization", action='store_true')
    parser.add_argument('-p', help="Run prediction", action='store_true')

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
                    " -noise "  + args.noise + \
                    " -i "      + infile 

        if args.out != None:
            ggp_arg +=  " -o " + args.out
            out_dir =  args.out
        else:
            out_dir = infile[:-4] + '_out'

        if args.m:
            ggp_arg += ' -m '
        if args.p:
            ggp_arg += ' -p '

        run_it = True
        # deal with parameter files
        if len(args.parameters)>0:
            if args.parameters[0] == "infer":
                parameter_files = get_parameter_files(infile, args.out)
            else:
                parameter_files = args.parameters

            if len(parameter_files) == 2:
                if len(args.free_up)>0:
                    free_up_parameters(parameter_files[1], args.free_up)
            else:
                run_it = False
                    
        suffix_param_files = []
        if len(args.suffix)>0:
            parameter_files = []
            for suffix in args.suffix:
                parameter_file = infile[:-4] + '_'+ suffix + '.txt'
                parameter_files.append(parameter_file)
                suffix_param_files.append(parameter_file)
            
        elif len(args.parameters)>0:
            parameter_files = args.parameters

        if len(parameter_file)==2:
            run_it=False
            print("Not enough parameter files found -> don't run")

        ggp_arg +=  " -b " + get_arg_list(parameter_files)
        
        
        # ============ run! ============ #
        if run_it:
            if args.rerun:
                run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"], args.verbose)
            else:
                # prediction_file = look_for_prediction_file(out_dir, infile)
                output_file = look_for_output_file(out_dir, infile, "iterations")
                output_file_p = look_for_output_file(out_dir, infile, "prediction")

                # prediction file does not exist
                if output_file == None:
                    print("The input file", infile, "has no output file yet -> RUN")
                    run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"], args.verbose)
                elif args.p and output_file_p==None:
                    print("The input file", infile, "has no prediction file yet -> RUN")
                    run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"], args.verbose)

                # prediction file is older than one of the paramfiles
                else:
                    if len(suffix_param_files)>0:
                        is_old = os.path.getmtime(output_file) < os.path.getmtime(suffix_param_files[0]) or os.path.getmtime(output_file) < os.path.getmtime(suffix_param_files[1])
                    else:
                        is_old = os.path.getmtime(output_file) < os.path.getmtime(parameter_files[0]) or os.path.getmtime(output_file) < os.path.getmtime(parameter_files[1])
                        
                    if is_old:
                        print(output_file, " is older than parameter files -> RUN")
                        run_command(config["bin"] + ' ' + ggp_arg, args.dryrun, config["iscluster"], args.verbose)
                    
                    # output file is up to date
                    else:
                        print(output_file, "is already there and up-to-date! -> No need to run")

# ================================================================================ #
if __name__ == "__main__":
    main()