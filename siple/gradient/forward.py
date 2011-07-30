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

import numpy
from siple.reporting import msg

class LinearForwardProblem:
  """
  We are interested in solving an ill-posed problem
  
  .. math::          T(x) = y
         
  where :math:`T:X\\rta Y` is a linear map between Hilbert spaces.  One class of methods
  for solving such a problem involve iteratively methods that approximately minimize
  
  .. math::  J(x)  = || y - T(x) ||^2_Y
        
  These techniques require the following maps: 
  
      The forward problem: T

      the adjoint of T:    T*

  This class encapsulates the solution of all these maps for a single forward problem,
  as well as defining the norms on the domain and range spaces.
  """

  def T(self,x,out=None):
    """
    Implements the forward linear problem.  Returns T(x) in the variable :data:`out`.
    """
    raise NotImplementedError()

  def TStar(self,y,out=None):
    """
    Implements the adjoint of the forward linear problem.  Returns T^*(y) in the variable :data:`out`.
    """
    raise NotImplementedError()

  def domainIP(self,u,v):
    """
    Implements inner product in the domain X.
    """
    raise NotImplementedError()
    
  def rangeIP(self,u,v):
    """
    Implements inner product in the range Y.
    """
    raise NotImplementedError()

  def testTStar(self, d, r):
    """
    Helper method for determining the adjoint has been coded correctly.

    Returns the following inner products:

      <T(d),r>_Y

      <d,T*(r)>_X

    If the linearization is correct and its adjoint is correct, then these should be the same number.
    """
    Td = self.T(d)
    TStarR = self.TStar(r)

    return (self.domainIP(d,TStarR), self.rangeIP(Td,r))


class NonlinearForwardProblem(LinearForwardProblem):
  """
  We are interested in solving an ill-posed problem
  
         F(x) = y
         
  where F:X->Y is a nonlinear map between Hilbert spaces.  One class of methods
  for solving such a problem involve iteratively methods that approximately minimize
  
    J(x)  = 1/2 || y - F(x) ||^2_Y
        
  These techniques require the following maps: 
  
      the forward problem: F
      its linearization:   T
      the adjoint of T:    T*
      
  This class encapsulates the solution of all these maps for a single forward problem,
  as well as defining the norms on the domain and range spaces.
  """
  
  def F(self, x,out=None,guess=None):
    """
    Returns the value of the forward problem at x.
    
    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.
    
    Storage in 'out', if given, is used for the return value.
    """
    raise NotImplementedError()
    
  def T(self,d,out=None):
    """
    Returns the value of the linearization, T, of F, at the point x specified previously in linearizeAt, 
    in the direction d.
    
    Storage in 'out', if given, is used for the return value.
    """
    raise NotImplementedError()

  def TStar(self,r,out=None):
    """
    Let T be the linearization of F at the point x (at the point x specified previously in linearizeAt).  
    Its adjoint is T*.  This method returns the value of T* in the direction r.
    
    Storage in 'out', if given, is used for the return value.
    """
    raise NotImplementedError()

  def linearizeAt(self,x,guess=None):
    """
    Instructs the class that subsequent calls to T and TStar will be conducted for the given value of x.

    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.
    """
    raise NotImplementedError()

  def evalFandLinearize(self,x,out=None,guess=None):
    """
    Computes the value of F(x) and locks in a linearization.  Sometimes there are efficiencies that
    can be acheived this way.
    
    Default implementation simply calls F, then linearizeAt.
    """
    Fx = self.F(x,out=out,guess=guess)
    self.linearizeAt(x,guess=Fx)
    return Fx

  def rangeIP(self,a,b):
    """
    Computes the inner product of two vectors in the range.
    """
    raise NotImplementedError()

  def domainIP(self,a,b):
    """
    Computes the inner product of two vectors in the domain.
    """
    raise NotImplementedError()

  def testT(self, x, h, t=1e-4,guess=None):
    """
    Helper method for determining the the forward linearziation has been coded correctly.
    
    Returns the difference quotient (F(x+th)-F(x))/t and the value T(h).  If the linearization
    has been coded correctly, these vectors should be close to each other.
    """
    xp = x.copy()
    xp.axpy(t,h)

    y=self.F(x,guess=guess)
    yp = self.F(xp,guess=y)

    Fp = yp.copy()
    Fp.axpy(-1,y)
    Fp /= t

    self.linearizeAt(x,guess=y)
    Fp2 = self.T(h)
    
    return (Fp,Fp2)

  def testTStar(self, x, d, r, guess=None):
    """
    Helper method for determining the adjoint has been coded correctly.
  
    Returns the following inner products:
    
      <T(d),r>_Y
      <d,T*(r)>_X
      
    If the linearization is correct and its adjoint is correct, then these should be the same number.
    """
    self.linearizeAt(x,guess=guess)

    Td = self.T(d)
    TStarR = self.TStar(r)

    return (self.domainIP(d,TStarR), self.rangeIP(Td,r))



class ForwardProblemLineSearchAdaptor:
  """
  Let F:X->Y be the nonlinear map between Hilbert spaces defined in a ForwardProblem. 
  
  We are interested in approximately minimizing a functional of the form
  
      J(x) = 0.5 * || y-F(x) ||_Y^2
  
  Doing this frequently involves a line-search in the direction 'd', i.e. the function
  
      Phi(t) = J(x+t*d)

  is minimized over t>=0.
  
  This class encapsulates the conversion of a ForwardProblem into the map Phi defined above.  
  """
  
  def __init__(self,forward_problem):
    """
    Initialization from a ForwardProblem
    """
    self.forward_problem = forward_problem
    self.z = None
    self.r = None
    self.Fz = None
    self.Td = None

  def eval(self,x,d,y,t):
    """
    Efficiently evaluate the function Phi(t) and Phi'(t) where
    
        Phi(t) = 0.5 * || y- F(x+td) ||_Y^2

    as described in the classlevel documentation.
    """

    if self.z is None:
      self.z = x.copy()
      self.r = y.copy()
      # self.Fz = dolfinutils.vectorLike(y)
      self.Td = y.vector_like()
    else:
      self.z.set(x)
      self.r.set(y)

    forward_problem = self.forward_problem

    self.z.axpy(t,d)

    try:
      self.Fz = forward_problem.evalFandLinearize(self.z,out=self.Fz,guess=self.Fz)
    except Exception as e:
      msg('Exception during evaluation of linesearch at t=%g;\n%s',t,str(e))
      return (numpy.nan,numpy.nan,None)
    self.Td=forward_problem.T(d,out=self.Td)

    self.r -= self.Fz
    J = 0.5*forward_problem.rangeIP(self.r,self.r)
    Jp = -self.forward_problem.rangeIP(self.Td,self.r)

    # Extra return value 'None' can be used to store extra data.
    return (J,Jp,self.Fz.copy())
