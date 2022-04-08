#!/bin/bash

#SBATCH --job-name=estim_param
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --time=23:00:00
#SBATCH --qos=1day

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

module load matplotlib
python3 estimate_init_parameters_theo.py -d ../../experimental_data/data_theo_input_20220325/experimental_data_GFP -s GM.1_parameters_GFP SPM.1_parameters_GFP