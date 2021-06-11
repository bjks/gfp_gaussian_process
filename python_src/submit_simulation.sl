#!/bin/bash

#SBATCH --job-name=cells_simulation
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --time=04:00:00
#SBATCH --qos=6hours

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=bjoern.kscheschinski@unibas.ch

python3 cells_simulation.py --o test --n 110 --l 20 --d 10 --i /scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210401_constitexpr_inference --g s
