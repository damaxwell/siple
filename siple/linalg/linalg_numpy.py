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

from linalg_abstract import AbstractVector
import numpy as np

class NumpyVector(AbstractVector):
  """Implements the siple.linalg.AbstractVector interface for a numpy array."""
  def __init__(self,u):
    if isinstance(u,tuple):
      u = np.ndarray(u)
    elif not isinstance(u,np.ndarray):
      raise ValueError("An NumpyVector can only be constructed from a numpy array or a size specification: found %s" % u)

    self._core = u

  def _set_from_abstract(self,rhs):
    self._core[:] = rhs._core[:]

  def _set_from_array(self,rhs):
    self._core[:] = rhs[:]

  def acc(self,rhs):
    self._core += rhs._core
  
  def scale(self,t):
    self._core *= t
  
  def axpy(self,t,v):
    self._core += t*v._core
  
  def copy(self):
    return NumpyVector(self._core.copy())
  
  def vector_like(self):
    return NumpyVector(np.ndarray(self._core.shape))
    
  def zero_like(self):
    z = np.zeros(self._core.shape)
    return NumpyVector(z)
  
  def dim(self):
    return self._core.shape[0]
  
  def core(self):
    return self._core

  def norm(self,name):
    if name == 'l2':
      return np.linalg.norm(self._core,2)
    if name == 'l1':
      return np.linalg.norm(self._core,1)
    if name == 'linf':
      return np.linalg.norm(self._core,np.inf)
  
  def __repr__(self):
    return self._core.__repr__()
  def __str__(self):
    return self._core.__str__()