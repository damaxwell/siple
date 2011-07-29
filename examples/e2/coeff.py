import siple
from siple.gradient.forward import NonlinearForwardProblem
from siple.linalg.linalg_numpy import NumpyVector
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from siple.reporting import pause, endpause
from matplotlib import pyplot as pp

class CoeffForwardProblem(NonlinearForwardProblem):
  """Implements the forward problem of finding a perioodic function :math:`u` 
  on :math:`[0,1]` solving 
  
  .. math:: -u'' + e^{\beta} u = 1.
  """
  
  def __init__(self,N):
    """Forward problem with *N* subintervals."""
    L = 1.
    self.h = L/N
    self.N = N

    self.u = np.ndarray((N,))
    self.rhs = np.ones((N,))

    self.A  = sparse
    self.B  = sparse
    self.AA = sparse
    self.C  = sparse

    
    self.scratch = np.ndarray((self.N+1,))

  def rangeIP(self,x,y):
    """:math:`L^2` inner product"""
    return <x,y>

  def domainIP(self,x,y):
    """:math:`L^2` inner product"""
    return <x,y>

  def evalF(self,beta,out=None):
    # build BB(v,w) = \int -e^{\beta} u v w
    # solve AA = rhs
    # return out
    pass
    
  def linearizeAt(self,beta):
    # solve for u
    # build C(v,w) = \int -e^{\beta} u v w
    # 
    pass
    
  def evalFAndLinearizeAt(self,beta,out=None):
    pass


  def T(self,x,out=None):
    """Computes the antiderivative of *x* with mean zero."""
    x = x.core()
    scratch = self.scratch

    if out is None:
      out = NumpyVector(np.ndarray(x.shape))

    scratch[0] = 0
    for k in range(self.N):
      scratch[k+1] = scratch[k] + self.h*x[k]

    y=out.core()
    for k in range(self.N):
      y[k] = 0.5*(scratch[k] + scratch[k+1])
    
    ybar = sum(y)/self.N
    y -= ybar

    return out

  def TStar(self,r,out=None):
    """Computes the adjoint of :math:`T`, in this case :math:`-T`."""
    out = self.T(r,out=out)
    out.scale(-1)
    return out
