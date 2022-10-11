# README run_all_dany.py
```
usage: run_all_dany.py [-h] -d DIR -c CSV_CONFIG [-o OUT] [-s SUFFIX]
                       [-space SPACE] [-t TOL] [--dryrun] [--local] [-m] [-p]
                       [-j]

Run gfp_gaussian for all files in dir (written for Dany's data sets)

optional arguments:
  -h, --help     show this help message and exit
  -d DIR         directory with input files
  -c CSV_CONFIG  Csv config file
  -o OUT         Output dir (None)
  -s SUFFIX      Suffix of parameter files ("_parameters.txt")
  -space SPACE   Space, log or linear (log)
  -t TOL         Tolerance of maximization (1e-15)
  --dryrun       Shows what will be done
  --local        Do not submit job, but run directly
  -m             Run maximization
  -p             Run prediction
  -j             Run joints
  ```