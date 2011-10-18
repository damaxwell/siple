import siple
from siple.gradient.forward import NonlinearForwardProblem
from siple.gradient.nonlinear import BasicInvertNLCG, BasicInvertIGN
from siple.linalg.linalg_numpy import NumpyVector
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from siple.reporting import pause, endpause
from math import sqrt, exp, pi
from matplotlib import pyplot as pp

class CoeffForwardProblem(NonlinearForwardProblem):
  """Implements the forward problem of finding a periodic function :math:`u` 
  on :math:`[0,1]` solving 
  
  .. math:: -u'' + e^{\beta} u = 1.
  """
  
  def __init__(self,N):
    """Setup finite element forward problem with *N* subintervals."""
    L = 1.
    self.h = L/N; h = self.h
    self.N = N

    self.x = np.arange(N)*h
    
    # 5th order gaussian quadrature points and weights
    self.quadpts = np.array([-sqrt(3./5.),0,sqrt(3./5.)])
    self.quadwts = np.array([5./9,8./9,5./9])

    # Data for finite-element basis functions
    l = (1.+self.quadpts.copy())/2.
    r = (1.-self.quadpts.copy())/2.
    self.b = [ l, r] # values of basis functions at quad points
    self.bprime = [-1./self.h, 1./self.h] # derivatives of basis functions

    # Matrix for \int \phi_i' \phi_j'
    A = sparse.lil_matrix((N,N))
    eps = 1./h
    for k in range(N):
      A[k, (k-1) % N] = -eps
      A[k, (k+1) % N] = -eps
      A[k, k] = 2*eps
    self.A = A.tocsr()

    # Matrix for \int \phi_i \phi_j
    B  = sparse.lil_matrix((N,N))
    for k in range(N):
      B[k,k] = 2*h/3.
      B[k,(k+1) % N] = h/6.
      B[k,(k-1) % N] = h/6.
    self.B = B.tocsr()

    #  Variable that will eventually hold the system matrix.
    self.system = None

    #  Variable that will eventually hold the low-order terms matrix
    self.low_order = None

    # Right-hand side for main PDE is always 1
    self.Frhs = self.B*np.ones((N,))

    # Storage for linearized operator right-hand side.
    self.Trhs = np.zeros((N,))

    # Storage for current value of beta, u when working with linearizations
    self.u = np.zeros((N,))
    self.beta = np.zeros((N,))

  def assemble_low_order(self,beta):
    """Assemble the matrix :math:`M[i,j]=\int e^\beta \phi_i \phi_j."""
    
    N = self.N
    self.low_order = sparse.lil_matrix((N,N))
    low_order = self.low_order
    quadpts = self.quadpts
    b = self.b; bprime = self.bprime
    l = b[0]; r=b[1]

    low_order *= 0.
    for k in xrange(N):
      for p in xrange(len(quadpts)):
        betap = beta[k]*l[p] + beta[(k+1) % N]*r[p] # value at quad pt.
        Jxw = self.quadwts[p]*self.h/2. # extra /2 since quad formulas are on [-1,1]
        for i in xrange(2):
          for j in xrange(2):
            I = (i+k) % N;  J = (j+k) % N; 
            low_order[I,J] += Jxw * exp(betap) * b[i][p] * b[j][p]

    # On the first go-round, we haven't yet converted to CSR format, so
    # do this now.
    if low_order.getformat() != 'csr':
      self.low_order = low_order.tocsr()

  def assemble_Trhs(self,beta,u,h):
    """Compute :math:`rhs[j] = \int e^\beta u \phi[j]` """
    N = self.N
    quadpts = self.quadpts
    b = self.b;
    l = b[0]; r=b[1]

    rhs = self.Trhs
    rhs[:] = 0
    
    for k in xrange(N):
      for p in xrange(len(quadpts)):
        Jxw = self.quadwts[p]*self.h/2
        kk = (k+1) % N
        betap = beta[k]*l[p] + beta[kk]*r[p] # value at quad pt.
        up = u[k]*l[p] + u[kk]*r[p] # value at quad pt.
        hp = h[k]*l[p] + h[kk]*r[p] # value at quad pt.
        rhs[k]  -= Jxw*exp(betap)*up*hp*b[0][p]
        rhs[kk] -= Jxw*exp(betap)*up*hp*b[1][p]

  def rangeIP(self,x,y):
    """:math:`L^2` inner product"""
    return np.dot(x.core(),self.B*y.core())

  def domainIP(self,x,y):
    """:math:`L^2` inner product"""
    return np.dot(x.core(),self.B*y.core())

  def F(self,beta,out=None,guess=None):

    self.assemble_low_order(beta.core())

    # scipy 'csr' matrix doesn't support +=, so be inefficient
    self.system = self.A + self.low_order
    
    u=spsolve(self.system,self.Frhs)

    return NumpyVector(u)
    # 
    # if not out is None:
    #   out.core()[:] = u[:]
    # else:
    #   out = NumpyVector(u)
    # return out
    
  def linearizeAt(self,beta,guess=None):
    self.evalFAndLinearize(beta)

  def evalFAndLinearize(self,beta,out=None,guess=None):
    out = self.F(beta,out=out)

    # Copy beta and u for use in later assembly.
    self.u[:] = out.core()[:]
    self.beta[:] = beta.core()[:]
    
    return out

  def T(self,h,out=None):
    h = h.core()

    if out is None:
      out = NumpyVector(h.shape)

    self.assemble_Trhs(self.beta,self.u,h)
    w = spsolve(self.system,self.Trhs)

    return NumpyVector(w)
    
  def TStar(self,g,out=None):
    g = g.core()
    
    if out is None:
      out = NumpyVector(g.shape)
    
    v = spsolve(self.system,self.B*g)
    self.assemble_Trhs(self.beta,self.u,v)
    w = spsolve(self.B,self.Trhs)

    return NumpyVector(w)


if __name__ == '__main__':
  from optparse import OptionParser

  usage =  """Usage: %prog [options]

Example: %prog -N 100 -n 0.1"""

  parser = OptionParser(usage=usage)
  parser.add_option("-N","--node_count",type='int',default=40,
                    help='number of subintervals')
  parser.add_option("-n","--l_infty_noise",type='float',default=0.001,
                    help='standard deviation of added noise')
  parser.add_option("-a","--algorithm",type='choice',choices=['sd','nlcg','ign'],default='nlcg',
                    help="algorithm to use [sd,nlcg,ign]: sd=steepest descent, nlcg=nonlinear conjugate gradient (default), ign=incomplete Gauss-Newton")
  parser.add_option("-t","--test_linearization",action='store_true',
                    help='test linearization and adjoint (and exit)')
  parser.add_option("-d","--discrepancy_fraction",type='float',default=1.0,metavar="D",
                    help='remove the fraction D of the actual error (D<1 to overfit and D>1 to underfit)')

  (options, args) = parser.parse_args()

  N = options.node_count
  Linf_error = options.l_infty_noise

  forward_problem = CoeffForwardProblem(N)
  x = forward_problem.x
  beta=NumpyVector((N,))
  beta.core()[:] = np.sin(x*2*pi)


  if options.test_linearization:    
    h = siple.rand.random_vector(beta,scale=1.)
    
    (Fp1,Fp2) =  forward_problem.testT(beta,h,t=1e-6)
    dF = Fp1.copy(); dF -= Fp2
    print 'Relative T error: %g' % (dF.norm('linf')/Fp1.norm('linf'))
    
    g = siple.rand.random_vector(beta,scale=1.)
    (ip1,ip2) = forward_problem.testTStar(beta,h,g)
    print 'Relative T* error: %g' % (abs(ip1-ip2)/abs(ip1))

  else:

    u_true = forward_problem.F(beta)
    u = u_true.copy()
    u += siple.rand.random_vector(u,scale=Linf_error)
    L2_error = Linf_error

    if options.algorithm == 'ign':
      Solver = BasicInvertIGN
    else:
      Solver = BasicInvertNLCG
      
    params = Solver.defaultParameters()
    params.ITER_MAX = 10000
    if options.algorithm == 'sd':
      params.steepest_descent = True
    solver = Solver(forward_problem,params=params)

    beta0 = u_true.zero_like()
    (betac,uc) = solver.solve(beta0,u,options.discrepancy_fraction*L2_error)

    pp.clf()
    pp.subplot(1,2,1)
    pp.plot(x,u.core(),x,uc.core())
    pp.legend(('Measured function', 'Computed function'))
    pp.subplot(1,2,2)
    pp.plot(x,beta.core(),x,betac.core())
    pp.legend(('True $\\beta$', 'Computed $\\beta$'))
    pp.draw()
    pp.show()
    
    endpause()
