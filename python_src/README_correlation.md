# README Correlation
This readme explains how correlation functions are calculated from the input data. 


## 1. Calculate joint posteriors
To calculate the joint posteriors the flag `-j` needs to be added to the inference run `./gfp_gaussian -j <other input>`. This will create the file `<prefix>_joints.csv` inlcuding the joint posteriors.


## 2. Calculate correlation functions
```
usage: correlation_from_joint.py [-h] -d DIR [-o OUTPUT_DIR]
                                 [-k KEY [KEY ...]] [-dt DT [DT ...]]
                                 [-n_data N_DATA] [-delimiter DELIMITER]

Correlation from joint matrix

optional arguments:
  -h, --help            show this help message and exit
  -d DIR                directory with input files OR joint file
  -o OUTPUT_DIR         directory for output
  -k KEY [KEY ...]      Keywords marking the files for a given dt (["acetate",
                        "glycerol", "glucose", "glucoseaa"])
  -dt DT [DT ...]       dt corresponding to keys provided by -k ([18.75, 6, 3,
                        1.5])
  -n_data N_DATA        Maximal number of data points that the maximial dt can
                        be appart (200)
  -delimiter DELIMITER  Delimiter in filename ("_")
```
The script `correlation_from_joint.py` takes the directory where the `<prefix>_joints.csv` files are located (`-d`). In the same derectory, the corresponding `<prefix>_prediciton.csv` file must be located. An ouput directory can be given (`-o`) otherwise the input  directory is used. This script will produce `<prefix>_correlations.npz` files. This step can be run as a slurm job with several cores.

The filenames of the data sets MUST include one of the keywords provided by `-k`, this will pick the dt provided by `-dt`, which should correspond to the time between measurements. The keyword must be seperated by the delimiter set by `-d`.


## 3. Plot correlation functions
The `<prefix>__correlations.npz` can be plotted using the python notebooks `plot_correlations_<...>_ipynb`.