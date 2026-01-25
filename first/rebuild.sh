#!/bin/sh

export PYTHONPATH=$PYTHONPATH:../shared

#fbase is the directory path to the daily OIv2 files
#Edit dates in the python scripts for desired epoch

time python3 old_first.py > alpha  # build traditional climatology
time python3 new_first.py > beta   # build epoch climatology -- epochYYYY.nc
time python3 maps.py      > gamma  # plot the trend and harmonics for epochal climo

# Residuals from 30 years
time python3 new_residuals2.py > delta # compute the residuals from the climatology

# compute residuals for a decade after the 30 years used for climatology
time python3 old_residuals.py  > epsi # traditional climatology 
time python3 new_res_decade.py > zeta # new climatology

# Recompute sumx1, sumx2, sumx3, sumx4 after subtracting off the traditional climatology 
#  -- same epoch. Numerics check.
time python3 trad_residuals1.py > eta

