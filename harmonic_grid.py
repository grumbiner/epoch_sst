from math import *
import numpy as np
import scipy
    
""" 
processing for simple time series of grids 
separate from processing for gridded data in a series due to matrix operations

""" 

# Summation functions -- cross terms for non-orthogonality of harmonics
def sssum(freq1, freq2, n):
  sum = 0. 
  for i in range(0, n):
    sum += sin(freq1*i)*sin(freq2*i)
  return sum
      
def scsum(freq1, freq2, n):
  sum = 0. 
  for i in range(0, n):
    sum += sin(freq1*i)*cos(freq2*i)
  return sum
  
def ccsum(freq1, freq2, n):
  sum = 0.
  for i in range(0, n):
    sum += cos(freq1*i)*cos(freq2*i)
  return sum

def sumsin(freq, t, n):
  sum = 0.
  for i in range(0,n):
    sum += sin(freq*t[i])
  return sum

def sumcos(freq, t, n):
  sum = 0.
  for i in range(0,n):
    sum += cos(freq*t[i])
  return sum

#***********************************************************----------!!
# x is the data vector, a, b are cos, sin amplitudes, respectively
# this changes in summations for grids

# time ranges 0 to t
def harmonic_coeffs(coeff, omega, t, m):
  for i in range(0, m):
    for j in range(0, m):
      coeff[2*i  , 2*j  ]  = ccsum(omega[i], omega[j], t)
      coeff[2*i+1, 2*j  ]  = scsum(omega[j], omega[i], t)
      coeff[2*i  , 2*j+1]  = scsum(omega[i], omega[j], t)
      coeff[2*i+1, 2*j+1]  = sssum(omega[i], omega[j], t)

def harmonic_solve(coeff, y, a, b, m):
# least squares solution for full generality:
  for i in range(0, y.shape[0]):
    for j in range(0, y.shape[1]):
      z = scipy.linalg.solve(coeff,y[i,j,:])
      for f in range(0,m):
        a[i,j,f] = z[2*f]
        b[i,j,f] = z[2*f+1]


