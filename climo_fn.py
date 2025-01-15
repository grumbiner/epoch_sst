import copy

def climo(intercept, slope, ampl, phase, freq, epoch, tag):
  delta = (tag - epoch).days
  sst = copy.deepcopy(intercept)
  sst += slope*delta
  for j in range(0,3):
    sst += ampl[j]*np.cos(phase[j] + freq[j]*delta)

  return sst

"""
Read in first pass -- 
  intercept, trend, harmonics 1-3 and their phase
fn to compute Tclim(tau) given the above (tau = days since 1 Sep 1981)

Second pass -- 
  Read in daily analyses
    subtract climatology
    accumulate stats on residuals
    accumulate terms for orthogonalizing w.r.t. Nino3.4
    write out residual field for the day
  Write out statistics on residuals
  Write out statistics on Nino3.4 orthogonalizing
  Maps of deviation norms
  Maps of correlations to Nino3.4

Third pass --
  Read in daily analyses
    subtract climatology
    accumulate stats on residuals
    write out residual field for the day
  Write out stats on residuals
  Maps of deviation norms

Offline:
  Compute + map:
    %variance explained by mean, trend, harmonics, Nino3.4
    Magnitude residual variance

"""
 
