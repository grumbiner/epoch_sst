#Edit dates for desired epoch

time python3 old_first.py > alpha  # build traditional climatology
time python3 new_first.py > beta   # build epoch climatology -- first_pass.nc
time python3 maps.py      > gamma  # plot the trend and harmonics for epochal climo
  
  # derive residual fields for traditional climatology

time python3 res1_nino34.py  # derive residual fields for epoch climatology
                             #   derive nino34 orthogonalization
 # plot the nino34 correlation and amplitude
