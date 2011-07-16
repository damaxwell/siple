import siple
from siple.gradient.forward import LinearForwardProblem
from siple.gradient.linear import  BasicKrylovCGNE
from siple.linalg.linalg_numpy import NumpyVector
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from siple.reporting import pause, endpause

class ForwardProblem(LinearForwardProblem):
  def __init__(self,L,N):
    if (N%2) == 0:
      raise ValueError('Discrete forward problem is singular for even values of \'N\'.\nUse an odd value of N, not %d' % N)
    h=float(L)/N
    A = sparse.lil_matrix((N+1,N+1))
    for k in range(N):
      A[k, (k+1) % N] = 0.5
      A[k, (k-1) % N] =  -0.5
      A[k,N] = -h
      A[N,k] = h
    A[N,N] = 1.
    self.A = A.tocsr()
    
    B  = sparse.lil_matrix((N,N))
    for k in range(N):
      B[k,k] = h/3.
      B[k,(k+1) % N] = h/6.
      B[k,(k-1) % N] = h/6.
    self.B = B.tocsr()

    self.rhs = np.zeros(N+1)

  def rangeIP(self,x,y):
    return np.dot(x.core(),self.B*y.core())
  
  def domainIP(self,x,y):
    return self.rangeIP(x,y)
  
  def T(self,x,out=None):
    # no reuse of memory is possible with numpy solve, so
    # we ignore 'out'.
    self.rhs[:-1] = self.B*x.core()
    self.rhs[-1] = 0
    y=spsolve(self.A,self.rhs)
    return NumpyVector(y[:-1])

  def TStar(self,r,out=None):
    x=self.T(r,out)
    x.scale(-1)
    return x

L=10
N=101;
h=float(L)/N
p = h*np.linspace(0,N-1,N)

y = NumpyVector(3*np.sin(p*np.pi*2/L)+np.cos(4*p*np.pi*2/L))
Linf_error = 0.01
L2_error = np.sqrt(L*Linf_error**2)

y += siple.rand.random_vector(y,scale=Linf_error)

x0 = y.zero_like()
params = BasicKrylovCGNE.defaultParameters()
params.steepest_descent = True
solver = BasicKrylovCGNE(ForwardProblem(L,N),params=params)
(xc,yc) = solver.solve(x0,y,L2_error)

from matplotlib import pyplot as pp
pp.plot(p,y.core(),p,yc.core(),p,xc.core())
pp.legend(('Measured function', 'Computed function', 'Computed derivative'))
pp.draw()

endpause()
