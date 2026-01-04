#!/bin/sh

export PYTHONPATH=$PYTHONPATH:../shared

#fbase is the directory path to the daily OIv2 files
#Edit dates in the python scripts for desired epoch

time python3 old_first.py > alpha  # build traditional climatology
time python3 new_first.py > beta   # build epoch climatology -- epochYYYY.nc
time python3 maps.py      > gamma  # plot the trend and harmonics for epochal climo

# compute residuals for a decade after the 30 years used for climatology
time python3 old_residuals.py > delta # traditional climatology 
time python3 new_res_decade.py > epsi # new climatology

 
