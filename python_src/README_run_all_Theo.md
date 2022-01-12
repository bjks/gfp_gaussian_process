# Readme for run_all_Theo script 
Run all data sets in specifc directory at once 

## Usage 
```
run_all_Theo.py [-h] -d DIR -b PARAMETERS [PARAMETERS ...] -c
                       CSV_CONFIG [-space SPACE] [-t TOL] [--dryrun] [--local]
                       [--fallback] [-m] [-p]
```

## Options 
```
optional arguments:
  -h, --help            show this help message and exit
  -d DIR                directory with inpput files
  -b PARAMETERS [PARAMETERS ...]
                        Parameter file(s)
  -c CSV_CONFIG         Csv config file
  -space SPACE          Space, log(default) or linear
  -t TOL                Tolerance of maximization
  --dryrun              Shows what will be done
  --local               Do not submit job, but run directly
  --fallback            Indicate that the parameter files are fallbacks if
                        there are no specific files
  --newparamfile        Creates third parameter file ('segment2') that is
                        identical to the segment1 but contains the beta from
                        segment0
  -m                    Run maximization
  -p                    Run prediction

```
The options `-c`, `-space`, `-t`, `-m`, `-p` behave as they do in direct `gfp_gaussian` runs. The `--dryrun` option shows what will be done without running anything. The other options are explained below.


## Input files
The option `-d` specifies the directory in which the input files are located. This will run all input files that end with `.csv`.

## Specifing parameter files
#### Default behaviour
In the default mode the parameter files that are submitted will be used for all input files. For example `-b parameters1.txt parameters2.txt` will use those parameter files for each job. 

#### Fallback option
If the option `--fallback` is set, those parameter files are only used for those data sets that do not have a specific one. That means, for each input file, for example `foo/bar.csv`, the script will first look for parameter files that are named `foo/bar_parameters1.txt` and `foo/bar_parameters2.txt` and use those instead of the fallback files `parameters1.txt` and `parameters2.txt`. The motivation of that behaviour is that, one can start running all data sets at using the default parameter files and then write specific ones, for cases where the initial parameters need to be modified that will be automatically preferred in the next run.

#### Infer parameter files
If the parameter file flag is ste to `infer`, each data set will be run using the generated parameter file that was generated by previous maximization runs. This is meant to be used for predictions. Additionally, the option `newparamfile` creates third parameter file (named `segment2`) that is identical to the `segment1` but contains the `beta` from `segment0`. This takes into account that the aquisition time (and thus the bleaching rate) does not change when the switch happens, but is changed later. 


## Running on the cluster (default) or locally
Unless the `--local` option is used, the script will (try to) submit jobs to slurm. For that it uses the template `submit_ggp_run.sl` which is in the same folder:
```
#!/bin/bash

#SBATCH --job-name=ggp
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --time=23:00:00
#SBATCH --qos=1day

#SBATCH --output=std.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

${COMMAND}
```
All runs will submitted as individual jobs. Of course the slurm script can be modified if needed.

If the script is run locally, everything is run one after another directly.