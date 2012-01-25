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

from siple import Parameters
from siple.reporting import msg
from siple.exceptions import IterationCountFailure
from math import sqrt

class KrylovSolver:
  """
  Base class for solving a linear ill-posed problem
    
    T(x) = y
    
  via an iterative Krylov method.
  """

  @staticmethod
  def defaultParameters():
    return Parameters('KrylovSolver', ITER_MAX=200, mu=1.1, cg_reset=0, steepest_descent=False)

  def __init__(self, params=None):
    self.params = self.defaultParameters()
    if not params is None: self.params.update(params)

    self.iteration_listeners = []

  def forwardProblem(self):
    """
    Returns the LinearForwardProblem that defines the inverse problem. 
    """
    raise NotImplementedError()

  def solve(self,x0,y,*args):
    """
    Run the iterative method starting from the initial point x0.
    """  
    raise NotImplementedError()

  def addIterationListener(self,listener):
    """
    Add an object to be called after each iteration.
    """
    self.iteration_listeners.append(listener)

  def iterationHook(self,count,x,y,d,r,*args):
    """
    Called during each iteration with the pertinent computations.  Handy for debugging and visualization.
    """
    for listener in self.iteration_listeners:
      listener(self,count,x,y,d,r,*args)

  def initialize(self,x0,y,*args):
    """
    This method is a hook called at the beginning of a run.  It gives an opportunity for the class to 
    set up information needed to decide conditions for the final stopping criterion.
    
    It may also be that the initial data 'x0' expresses the the initial data for the problem T(x)=y
    indirectly. Or it could be that x0 and y are expressed as some sort of concrete vectors rather than
    some invtools.AbstractVector's.
    So initialize returns a pair of AbstractVectors (x0,y) which are possibly modified and or wrapped
    versions of the input data.
    
    The arguments \*args are passed directly from 'run'.
    """
    return (x0,y)

  def finalize(self,x,y):
    """
    This method is a hook called at the end of a run, and gives the class a way to make adjustments to x and y before
    finishing the run.
    """
    return (x,y)

  def stopConditionMet(self,iter,x,y,r):
    """
    Given a current iteration number, current value of x, desired value y of F(X), and current residual, 
    returns whether the stop condition has been met.
    """
    raise NotImplementedError()


class KrylovCG(KrylovSolver):
  """
  Base class for solving an ill-posed linear problem
    
    T(x) = y
    
  using a conjugate gradient method.  Necessarily, T:X->X must be self adjoint.
  """


  def solve(self,x0,y,*args):
    (x,y) = self.initialize(x0,y,*args)

    cg_reset = x.size()
    if( self.params.cg_reset != 0): cg_reset = self.params.cg_reset

    forward_problem = self.forwardProblem()
    r=y.copy()
    Tx = forward_problem.T(x)
    r -= Tx

    normsq_r = forward_problem.domainIP(r,r)

    # d = r
    d = r.copy()

    # Eventually will contain T(d)
    Td = None

    count = 0
    while True:
      if count > self.params.ITER_MAX:
        raise IterationCountFailure(self.params.ITER_MAX)
      count += 1

      if self.stopConditionMet(count,x,y,r):
        msg('done at iteration %d', count)
        break

      if self.params.verbose:
        msg('solving linear problem')
      Td = forward_problem.T(d,out=Td)

      self.iterationHook( count, x, y, d, r, Td )

      alpha = normsq_r/forward_problem.domainIP(d,Td)
      if( (count+1 % cg_reset) == 0 ): 
         msg('resetting cg via steepest descent')
         alpha = 1
      
      # x = x + alpha*d
      x.axpy(alpha,d)

      # r = r - alpha*Td
      r.axpy(-alpha,Td)

      prev_normsq_r = normsq_r
      normsq_r = forward_problem.domainIP(r,r)
      beta = normsq_r / prev_normsq_r
      if( (count+1 % cg_reset) == 0 ): beta = 0
      
      if(self.params.steepest_descent):
        beta = 0
      # d = T*r + beta*d
      d *= beta
      d += r


    y = forward_problem.T(x)
    return self.finalize(x,y)
      
      
class KrylovCGNE(KrylovSolver):  
  """
  Base class for solving an ill-posed linear problem
    
    T(x) = y
    
  using a conjugate gradient method applied to the normal equation

    T^*T x = T^* y
  """

  def solve(self,x0,y,*args):
    (x,y) = self.initialize(x0,y,*args)

    forward_problem = self.forwardProblem()    
    Tx = forward_problem.T(x)
    r = y.copy()
    r -= Tx

    cg_reset = x.size()
    if( self.params.cg_reset != 0): cg_rest = self.params.cg_reset

    TStarR = forward_problem.TStar(r)
    normsq_TStarR = forward_problem.domainIP(TStarR,TStarR)

    # d = T^* r
    d = TStarR.copy()

    # Eventual storage for T(d)
    Td = None

    count = 0
    while True:
      if count > self.params.ITER_MAX:
        raise IterationCountFailure(self.params.ITER_MAX)
        break
      count += 1

      if self.stopConditionMet(count,x,y,r):
        msg('done at iteration %d', count)
        break

      Td = forward_problem.T(d,out=Td)

      self.iterationHook( count, x, y, d, r, Td, TStarR )

      alpha = normsq_TStarR/forward_problem.rangeIP(Td,Td)
      if( (count+1 % cg_reset) == 0 ): 
        msg('resetting cg via steepest descent')
        alpha = 1

      # x = x + alpha*d
      x.axpy(alpha,d)

      # r = r - alpha*Td
      r.axpy(-alpha,Td)

      # beta = ||r_{k+1}||^2 / ||r_k||^2
      prev_normsq_TStarR = normsq_TStarR
      TStarR = forward_problem.TStar(r,out=TStarR)
      normsq_TStarR = forward_problem.domainIP(TStarR,TStarR)
      beta = normsq_TStarR/prev_normsq_TStarR
   
      if( (count+1 % cg_reset) == 0 ): beta = 0

      if(self.params.steepest_descent):
        beta = 0

      # d = T*r + beta*d
      d *= beta
      d += TStarR


    Tx = forward_problem.T(x, out=Tx)
    return self.finalize(x,Tx)

class BasicKrylovCG(KrylovCG):
  """
  Implements the CG regularization method for solving the linear ill posed problem
  
    T(x) = y
    
  using the Morozov discrepancy principle.  The discrepancy of 'x' is
  
    ||y-T(x)||_Y
  
  and the algorithm is run until a target discrepancy (specified as an argument to solve)
  is reached.
  
  The specific problem to solve is specified as an argument to the constructor.
  """
  def __init__(self,forward_problem,params=None):
    KrylovCG.__init__(self,params)
    self.forward_problem = forward_problem
    
  def forwardProblem(self):
    """
    Returns the LinearForwardProblem that defines the inverse problem. 
    """
    return self.forward_problem

  def solve(self,x0,y,targetDiscrepancy):
    """
    Run the iterative method starting from the initial point x0.
    
    The third argument is the desired value of ||y-T(x)||_Y
    """
    return KrylovCG.solve(self,x0,y,targetDiscrepancy)

  def initialize(self,x0,y,targetDiscrepancy):
    """
    This method is a hook called at the beginning of a run.  It gives an opportunity for the class to 
    set up information needed to decide conditions for the final stopping criterion.

    It may be that the initial data 'x0' expresses the the initial data for the problem T(x)=y
    indirectly. Or it could be that x0 and y are expressed as dolfin.Function's rather than dolfind.GenericVectors.
    So initialize returns a triple of vectors (x0,y) which are possibly modified versions of the input data.

    The arguments \*args are passed directly from 'run'.
    """
    self.targetDiscrepancy = targetDiscrepancy
    return (x0,y)

  def stopConditionMet(self,iter,x,y,r):
    """
    Given a current iteration number, current value of x, desired value y of F(X), and current residual, 
    returns whether the stop condition has been met.
    """
    return sqrt(self.forward_problem.rangeIP(r,r)) <= self.params.mu*self.targetDiscrepancy

class BasicKrylovCGNE(KrylovCGNE):
  """
  Implements the CGNE regularization method for solving the linear ill posed problem

    T(x) = y

  using the Morozov discrepancy principle.  The discrepancy of 'x' is

    ||y-T(x)||_Y

  and the algorithm is run until a target discrepancy (specified as an argument to solve)
  is reached.

  The specific problem to solve is specified as an argument to the constructor.
  """
  def __init__(self,forward_problem,params=None):
    KrylovCGNE.__init__(self,params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    """
    Returns the LinearForwardProblem that defines the inverse problem. 
    """
    return self.forward_problem

  def solve(self,x0,y,targetDiscrepancy):
    """
    Run the iterative method starting from the initial point x0.

    The third argument is the desired value of ||y-T(x)||_Y
    """
    return KrylovCGNE.solve(self,x0,y,targetDiscrepancy)

  def initialize(self,x0,y,targetDiscrepancy):
    """
    This method is a hook called at the beginning of a run.  It gives an opportunity for the class to 
    set up information needed to decide conditions for the final stopping criterion.

    It may be that the initial data 'x0' expresses the the initial data for the problem T(x)=y
    indirectly. Or it could be that x0 and y are expressed as dolfin.Function's rather than dolfind.GenericVectors.
    So initialize returns a triple of vectors (x0,y) which are possibly modified versions of the input data.

    The arguments \*args are passed directly from 'run'.
    """
    self.targetDiscrepancy = targetDiscrepancy
    return (x0,y)

  def stopConditionMet(self,iter,x,y,r):
    """
    Given a current iteration number, current value of x, desired value y of F(X), and current residual, 
    returns whether the stop condition has been met.
    """
    disc = sqrt(self.forward_problem.rangeIP(r,r))
    target = self.params.mu*self.targetDiscrepancy
    if self.params.verbose:
      msg('Iteration %d: discrepancy %g target %g',iter,disc,target)
    return disc <= target
