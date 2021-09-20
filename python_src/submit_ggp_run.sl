#!/bin/bash

#SBATCH --job-name=ggp
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --time=23:00:00
#SBATCH --qos=1day

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

${COMMAND}