import siple
from siple.gradient.forward import LinearForwardProblem
from siple.gradient.linear import  BasicKrylovCG
from siple.linalg.linalg_numpy import NumpyVector
from siple.reporting import pause
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from matplotlib import pyplot as pp

class ForwardProblem(LinearForwardProblem):
  def __init__(self,L,N):
    h=float(L)/N

    A = sparse.lil_matrix((N+1,N+1))
    eps = 1/h/h
    for k in range(N):
      A[k, (k-1) % N] = -eps
      A[k, (k+1) % N] = -eps
      A[k, k] = 2*eps
      A[k,N] = h
      A[N,k] = h
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
    return x

class Solver(BasicKrylovCG):
  def iterationHook(self,k,x,y,d,r,*args):
    pp.plot(x.core())
    pp.draw()

L=10
N=100;
h=float(L)/N
p = h*np.linspace(0,N-1,N)

y = NumpyVector(3*np.sin(p*np.pi*2/L)+np.cos(4*p*np.pi*2/L))
Linf_error = 0.1
L2_error = np.sqrt(L*Linf_error**2)

r=siple.rand.random_direction(y)
d=siple.rand.random_direction(y)
print ForwardProblem(L,N).testTStar(d,r)

y += siple.rand.random_vector(y,scale=Linf_error)

x0 = y.zero_like()
params = Solver.defaultParameters()
params.steepest_descent = True
solver = Solver(ForwardProblem(L,N),params)
(xc,yc) = solver.solve(x0,y,L2_error)

pp.clf()
pp.plot(p,y.core(),p,yc.core())
pp.show()