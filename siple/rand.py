############################################################################
#
#  This file is a part of siple.
#
#  Copyright 2010 David Maxwell
#
#  siple is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
# 
############################################################################

import numpy as np

# This module contains utility functions for generating vectors
# with randomly generated coefficients.

def random_vector(u,scale=1,seed=None):
  """Returns a vector in the same function space as u with normally
  distributed coefficients with mean zero and standard deviation 'scale'
  """
  if not (seed is None):
    np.random.seed(seed)

  h = u.vector_like()
  h.set(np.random.normal(scale=scale, size=u.size()))
  return h
  

def random_direction(u,scale=1,seed=None):
  """
  Returns a random vector with normally distributed entries with standard 
  deviation ||u||_\infty * scale.
  
  Useful for obtaining a random direction vector with size on the order of u.
  """
  scale *= u.norm('linf')
  
  return random_vector(u,scale=scale,seed=seed)


def addnoise(u,relative_percent_error,seed=None):
  """
  Add gaussian noise to a vector.  The standard deviation of the noise
  is (relative_percent_error)\% of the maximum value of the vector u.
  """
  if not (seed is None):
    np.random.seed(seed)

  noisy_u = u.copy()    
  max_u = u.norm('linf')

  scale = (relative_percent_error/100.*max_u)
  noisy_u += random_vector(u,scale=scale,seed=seed)

  return (nu,scale)

def noiseForGrid(N,M,dx,dy,minWaveLength):
  """
  Returns noise on an N by M grid with spacing dx and dy with a frequency spectrum as follows:

    1) frequencies corresponding to wavelengths small than minWaveLength are all zero
    2) amplitude of each nonzero component in normally distributed
    3) phase of each nonzero component is uniformly distributed

    The noise has zero mean and is scaled to have RMS=1.
  """
  x=np.arange(0,N,1)
  y=np.arange(0,M,1)

  (Y,X)=np.meshgrid(y,x)

  r = np.random.randn(N,M)

  theta = np.random.rand(N,M)*2*np.pi*(1J)
  freqs = r*np.exp(theta)

  a=N*dx/minWaveLength
  b=M*dy/minWaveLength

  mask1 = np.sqrt(X**2/a**2+Y**2/b**2)>1
  mask2 = np.sqrt((X-N)**2/a**2+Y**2/b**2)>1

  freqs[mask1&mask2]=0

  # This could be done faster with fancy slicing
  for i in range(1,N/2):
    for j in range(1,M/2):
      freqs[N-i,M-j] = np.conj(freqs[i,j])
      freqs[i,M-j] = np.conj(freqs[N-i,j])
  for i in range(1,N/2):
    freqs[N-i,0]=np.conj(freqs[i,0])
  for j in range(1,M/2):
    freqs[0,M-j]=np.conj(freqs[0,j])
  # Ensure we have mean zero
  freqs[0,0]=0

  noise = np.real(np.fft.ifft2(freqs))

  rms = np.sqrt(sum(sum(noise*noise))/(N*M))
  noise/=rms

  return noise
