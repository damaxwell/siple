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

class AbstractVector:
  """
  Defines the interface that a vector is required to implement for use in the
  inverse problem routines.  Typically a subclass will have a member variable
  that contains a class that is the concrete representation of the vector (e.g. 
  a FEniCS Function or a PETSc Vec) and translates the AbstractVector methods
  to operations on its concrete data.
  
  An abstract vector does not have a notion of of an assoiciated basis, so
  it does not make sense to access the 'coefficients' of an abstract vector.
  Hence there is no operator[] associated with an AbstractVector.
  
  When the underlying vector is needed, use the core() method to access it.
  """
  
  def set(self,rhs):
    """Sets the value of the vector to the value of 'rhs'. If 'rhs' is not of the same
    class as the current vector, an exception may be thrown.
    """
    if isinstance(rhs,AbstractVector):
      self._set_from_abstract(rhs)
    else:
      self._set_from_array(rhs)

  def _set_from_abstract(self,rhs):
    raise NotImplementedError()

  def _set_from_array(self,rhs):
    """Sets the vector from anything array-like. The interpretation of this
    operation is implementation dependent."""
    raise NotImplementedError()


  def acc(self,rhs):
    """Adds the vector 'rhs' to the current vector (i.e. accumulates rhs)."""
    raise NotImplementedError()

  def sum(self,rhs,storage=None):
    """Returns a vector that is the sum of the current vector and 'rhs'.  Uses 'storage'
    for the return value if it is specified."""
    if storage is None:
      storage = self.copy()
    else:
      storage.set(self)
    storage.acc(rhs)
    return storage

  def diff(self,rhs,storage=None):
    """Returns a vector that is the difference of the current vector and 'rhs'.  Uses 'storage'
    for the return value if it is specified."""
    if storage is None:
      storage = self.copy()
    else:
      storage.set(self)
    storage.axpy(-1,rhs)
    return storage

  def mul(self,t,storage=None):
    """Returns a vector that is the current vector, using 'storage' for the return value if given."""
    if storage is None:
      storage = self.copy()
    else:
      storage.set(self)
    storage.scale(t)
    return storage

  def lincomb(self,t0,t1,v,storage=None):
    """Returns the vector t0*u+t1*v where u is the current vector.  Uses 'storage'
    for the return value if it is specified"""
    if storage is None:
      storage = self.copy()
    else:
      storage.set(self)
    storage.scale(t0)
    storage.axpy(t1,v)
    return storage

  def dim(self):
    """Returns the dimension of the ambient vector space"""
    raise NotImplementedError()

  def size(self):
    return self.dim()

  def core(self):
    """Returns the underlying concrete vector."""
    raise NotImplementedError()

  def scale(self,t):
    """Multiplies the current vector by the scalar 't'."""
    raise NotImplementedError()

  def axpy(self,t,v):
    """Adds t*v to the current vector."""
    raise NotImplementedError()

  def copy(self):
    """Returns a new vector with the same contents as sthe current one."""
    raise NotImplementedError()
    
  def vector_like(self):
    """Returns a new vector in the same vector space as the current vector.  The value
    of the vector is not specified."""
    raise NotImplementedError()
    
  def zero_like(self):
    """Returns a zero vector in the same vector space as the current vector."""
    raise NotImplementedError()

  def norm(self,name):
    """Returns a norm of the vector.  The string name describes the desired
    norm (e.g. 'l1', 'l2', 'linf').  There is no guarantee that an AbstractVector
    has a notion of any particular norm. Use with caution."""
    raise NotImplementedError()

  def __imul__(self,t):
    self.scale(t)
    return self

  def __iadd__(self,v):
    self.acc(v)
    return self

  def __isub__(self,v):
    self.axpy(-1,v)
    return self
  
  def __idiv__(self,t):
    self.scale(1./t)
    return self
    
  def __sub__(self,other):
    return self.diff(other)
    
  def __add__(self,other):
    return self.sum(other)