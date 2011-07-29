import siple
from siple.gradient.forward import LinearForwardProblem
from siple.gradient.linear import  BasicKrylovCGNE
from siple.linalg.linalg_numpy import NumpyVector
import numpy as np
from siple.reporting import pause, endpause
from matplotlib import pyplot as pp

class AntiderivativeForwardProblem(LinearForwardProblem):
  """Implements the forward problem of finding the zero-mean antiderivative
  of a zero-mean function.
  """
  
  def __init__(self,L,N):
    """Forward problem with an interval of length *L* and *N* subintervals."""
    self.h = float(L)/N
    self.N = N
    self.L = L

    # Coordinates 
    self.p = self.h*np.linspace(0,N-1,N)+self.h/2

    self.scratch = np.ndarray((self.N+1,))
    
  def rangeIP(self,x,y):
    """:math:`L^2` inner product"""
    return self.h*np.dot(x.core(),y.core())
  
  def domainIP(self,x,y):
    """:math:`L^2` inner product"""
    return self.rangeIP(x,y)

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

class Solver(BasicKrylovCGNE):
  def iterationHook(self,count,x,y,d,r,*args):
    # p = self.forwardProblem().p
    # pp.clf()
    # pp.plot(d.core())
    # pp.draw()
    # pause()
    pass

  def solve(self,x,y,targetDisc):
    y = y.copy()
    ybar = sum(y)/len(y)
    y -= ybar
    
    x = x.copy()
    xbar = sum(x)/len(x)
    x -= xbar

    (xc,yc) = BasicKrylovCGNE.solve(self,NumpyVector(x),NumpyVector(y),targetDisc)

    xc.core()[:] += xbar
    yc.core()[:] += ybar
    
    return (xc.core(),yc.core())


if __name__ == '__main__':
  from optparse import OptionParser
  
  usage =  """Usage: %prog [options]

Example: %prog -L 10 -N 100 -n 0.1"""

  parser = OptionParser(usage=usage)
  parser.add_option("-L","--length",type='float',default=10,
                    help='interval length')
  parser.add_option("-N","--node_count",type='int',default=100,
                    help='number of subintervals')
  parser.add_option("-n","--l_infty_noise",type='float',default=0.1,
                    help='standard deviation of added noise')
  parser.add_option("-s","--steepest_descent",action='store_true',
                    help='use steepest descent rather than conjugate gradient')
  parser.add_option("-a","--test_adjoint",action='store_true',
                    help='test adjoint and exit')
  parser.add_option("-d","--discrepancy_fraction",type='float',default=1.0,metavar="D",
                    help='remove the fraction D of the actual error (D<1 to overfit and D>1 to underfit)')

  (options, args) = parser.parse_args()

  L = options.length
  N = options.node_count
  Linf_error = options.l_infty_noise

  forward_problem = AntiderivativeForwardProblem(L,N)

  h=float(L)/N

  if options.test_adjoint:
    y = NumpyVector((N,))
    d=siple.rand.random_vector(y,scale=1)
    d.core()[:] -= sum(d.core())/N

    r=siple.rand.random_vector(y,scale=1)
    r.core()[:] -= sum(r.core())/N

    (ip1,ip2) = forward_problem.testTStar(d,r)
    print 'Adjoint test:\n%f %f\nrelative error %f' % (ip1,ip1,abs(ip1-ip2)/abs(ip1))

  else:
    p = forward_problem.p
    y = 3*np.sin(p*np.pi*2/L)+np.cos(4*p*np.pi*2/L)+2.

    L2_error = np.sqrt(L*Linf_error**2)

    y += siple.rand.random_vector(NumpyVector(y),scale=Linf_error).core()

    params = Solver.defaultParameters()
    params.ITER_MAX = 10000
    params.steepest_descent = options.steepest_descent

    x0 = np.zeros(y.shape)
    solver = Solver(forward_problem,params=params)
    ( xc,yc) = solver.solve(x0,y,options.discrepancy_fraction*L2_error)

    from matplotlib import pyplot as pp
    pp.clf()
    pp.plot(p,y,p,yc,p,xc)
    pp.legend(('Measured function', 'Computed function', 'Computed derivative'))
    pp.draw()

    endpause()
