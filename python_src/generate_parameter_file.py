import argparse

def main():
    parser = argparse.ArgumentParser(
        description='Generate new parameter file from iteration file OR and final file of a gfp_gaussian run')

    parser.add_argument('--i',
                        dest='infile',
                        help='input file, iteration file from gfp_gaussian run',
                        required=True)
    parser.add_argument('--o',
                        dest="outfile",
                        help="output file, default: parameters_generated.txt",
                        required=False,
                        default="parameters_generated.txt")
    parser.add_argument('--s',
                        dest="step",
                        help="relative step size",
                        required=False,
                        default=0.01,
                        type=float)

    args = parser.parse_args()

    parameter_names =  ["mean_lambda", 
                        "gamma_lambda",
                        "var_lambda",
                        "mean_q",
                        "gamma_q",
                        "var_q",
                        "beta",
                        "var_x",
                        "var_g",
                        "var_dx",
                        "var_dg"]



    if args.infile.endswith("iterations.csv"):
        parameter_init = get_last_parameter_set(args.infile, parameter_names)
    elif args.infile.endswith("final.csv"):
        parameter_init = get_final_parameter_set(args.infile, parameter_names)
    else:
        print("File type (iterations/final) not recognized!\nQuit")
        return

    parameter_setting = get_parameter_setting(args.infile, parameter_names)

    print("Got parameter setting")
    print("----------------------------")
    for k,v in parameter_setting.items():
        print(f'{k : <13}', *v)


    print("\nNew init values")
    print("----------------------------")
    for k,v in parameter_init.items():
        print(f'{k : <13}', v)   

    write_param_file(parameter_setting, parameter_init, args.outfile, args.step)
    return

# ================================================================================ # 
# ================================================================================ #
def get_parameter_setting(filename, parameter_names):
    parameter_setting = dict.fromkeys(parameter_names)
    with open(filename) as fin:
        next(fin)
        for k in parameter_setting.keys():
            line = fin.readline()
            p_type = line.split(',')[2]
            bounds = line.split(',')[5:7]
            parameter_setting[k] = [p_type] + bounds
    return parameter_setting


def get_final_parameter_set(filename, parameter_names):
    final_set = []
    with open(filename) as fin:
        next(fin)
        for k in parameter_names:
            line = fin.readline()
            final_set.append(line.split(',')[-1].strip())
    print(final_set)
    return dict(zip(parameter_names, final_set))


def get_last_parameter_set(filename, parameter_names):
    current_set = []
    with open(filename) as fin:
        for _ in range(14):
            next(fin)
        while True:
            line = fin.readline()
            if not line:
                break
            try:
                float(line.split(',')[-1])
                current_set = line.split(',')[1:-1]
            except ValueError:
                pass
    return dict(zip(parameter_names, current_set))


# ================================================================================ #
def write_param_file(parameter_setting, parameter_init, outfile, step_scaling):
    with open(outfile, "w") as fin:
        fin.write("# Generated config file for simulated data\n")
        for k in parameter_setting.keys():

            step = str(max([float(parameter_init[k])*step_scaling , 1e-9]))
            if parameter_setting[k][0] == 'free':
                fin.write("{:s} = {:s}, {:s} \n".format(k, parameter_init[k], step) )
            elif parameter_setting[k][0] == 'fixed':
                fin.write("{:s} = {:s} \n".format(k, parameter_init[k]))
            else:
                fin.write("{:s} = {:s}, {:s}, {:s}, {:s} \n".format(k, 
                                                                parameter_init[k],
                                                                step,
                                                                parameter_setting[k][1],
                                                                parameter_setting[k][2]))


# ================================================================================ #
if __name__ == "__main__":
    main()