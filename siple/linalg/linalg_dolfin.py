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
import dolfin

class DolfinGenericVector(AbstractVector):
  """Implements the siple.linalg.AbstractVector interface for a dolfin generic vector."""
  def __init__(self,u):
    if isinstance(u,dolfin.Function):
      self._core = u.vector()
    elif isinstance(u,dolfin.GenericVector):
      self._core = u
    else:
      raise ValueError("An DolfinGenericVector only be constructed from a dolfin Function or a dolfin GenericVector: found a %s" % u)
  
  def _set_from_abstract(self,rhs):
    self._core[:] = rhs._core[:]

  def _set_from_array(self,rhs):
    self._core[:] = rhs[:]
  
  def acc(self,rhs):
    self._core += rhs._core
  
  def scale(self,t):
    self._core *= t
  
  def axpy(self,t,v):
    self._core.axpy(t,v.core())
  
  def copy(self):
    return DolfinGenericVector(self._core.copy())
  
  def vector_like(self):
    w =  self._core.factory().create_vector();
    w.resize(self._core.size())
    return DolfinGenericVector(w)
    
  def zero_like(self):
    z = self.vector_like()
    z._core.zero()
    return z
  
  def dim(self):
    return self._core.size()
  
  def core(self):
    return self._core

  def norm(self,name):
    return self._core.norm(name)

class DolfinFunctionVector(AbstractVector):
  """Implements the siple.linalg.AbstractVector interface for a dolfin Function."""

  def __init__(self,u):
    if isinstance(u,dolfin.Function):
      self._core = u
    elif isinstance(u,dolfin.FunctionSpace):
      self._core = dolfin.Function(u)
    else:
      raise ValueError("An DolfinFunctionSpaceVector wraps a dolfin Function or FunctionSpace: found a %s" % u)

  def set(self,rhs):
    self._core.vector()[:] = rhs._core.vector()[:]

  def acc(self,rhs):
    self._core.vector()[:] += rhs._core.vector()[:]

  def scale(self,t):
    self._core.vector()[:] *= t

  def axpy(self,t,v):
    self._core.vector().axpy(t,v.core().vector())

  def copy(self):
    rv = dolfin.Function(self._core.function_space())
    rv.vector()[:] = self._core.vector()[:]
    return DolfinFunctionSpaceVector(rv)

  def vector_like(self):
    w = dolfin.Function(self._core.function_space())
    return DolfinFunctionSpaceVector(w)

  def zero_like(self):
    z = self.vector_like()
    z._core.vector().zero()
    return z

  def dim(self):
    return self._core.vector().size()

  def core(self):
    return self._core

  def norm(self,name):
    return self._core.vector().norm(name)



def apply_form(F,x,y,scratch=None):
  vx = x
  if( isinstance(x,dolfin.Function) ):
    vx = x.vector();

  vy = y
  if( isinstance(y,dolfin.Function) ):
    vy = y.vector();

  if scratch is None:
    Fy = F*vy
  else:
    F.mult(vy,scratch)
    Fy = scratch
  return vx.inner(Fy)

def lhsVectorFor(A,out=None):
  outv=None
  if not out is None:
    outv = out.core()
  if (outv is None) or (outv.size() != A.size(0)):
    outv = A.factory().create_vector()
    outv.resize(A.size(0))
    out = DolfinGenericVector(outv)
  return out
