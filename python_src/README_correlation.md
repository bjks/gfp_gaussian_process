# README Correlation
This readme explains how correlation functions are calculated from the input data. 

The filenames of the data sets MUST include one of the following strings
```
"acetate_", "glycerol_",  "glucose_", "glucoseaa_"
```
This will set the time spacing in the second step.

## 1. Calculate joint posteriors
To calculate the joint posteriors the flag `-j` needs to be added to the inference run `./gfp_gaussian -j <other input>`. This will create the file `<prefix>_joints.csv` inlcuding the joint posteriors.


## 2. Calculate correlation functions
The filenames of the data sets MUST include one of the following strings
```
"acetate_", "glycerol_",  "glucose_", "glucoseaa_"
```
This will set the time spacing with is hardcoded at the moment (variable is called `dts` in `process_file()`) and should be identical to the measurment time spacing. 

The script `correlation_from_joint.py` takes the directory where the `<prefix>_joints.csv` files are located (`-d`). An ouput directory can be given (`-o`) otherwise the input  directory is used. This script will produce `<prefix>__correlations.npz` files. This step can be run as a slurm job with several cores.

## 3. Plot correlation functions
The `<prefix>__correlations.npz` can be plotted using the python notebooks `plot_correlations_<>_ipynb`.