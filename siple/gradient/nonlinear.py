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

import numpy as np
from siple.reporting import msg
from siple.exceptions import IterationCountFailure, NumericalFailure
from siple import Parameters
from siple.opt import linesearchHZ
from forward import ForwardProblemLineSearchAdaptor
from linear import KrylovSolver, BasicKrylovCGNE
from math import sqrt

class InvertIGN:
  """
  Implements solving an ill-posed minimization problem of the form
  
  .. math::       J(x) = B(y-\calF(x),y-\calF(x))

  where :math:`B` is a symmetric positive semidefinite bilinear form and :math:`\calF` is a (possibly nonlinear) function mapping between
  two Hilbert spaces.

  Minimization is done applying an incomplete Gauss-Newton algorithm until a stopping criterion is met.
  """

  @staticmethod
  def defaultParameters():
    """Default parameters that can subsequently be modified and passed to the constructor.

    :param ITER_MAX: maximum number of minimization iterations before failure
    :param mu: Morozov discrepancy principle scaling paramter
    :param cg_reset: number of iterations before a conjugate gradient reset (0=default)
    :param thetaMin: smallest acceptable goal improvment for linearized step
    :param thetaMax: maximum attempted goal improvment for linearized step
    :param kappaTrust: factor to shrink :data:`theta` by when linearized step fails
    :param rhoLow: linearized step fails when goal improves by less than :data:`rhoLow`*:data:`theta`
    :param rhoHigh: linearized step is very successful when goal improves by more than :data:`rhoHigh`*:data:`theta`
    :param verbose: True if extra logging output is desired
    """
    params = Parameters('InvIGN', ITER_MAX=200, mu=1.1, cg_reset=0, thetaMin=2**(-20), thetaMax=0.5, 
                                        kappaTrust=0.5, rhoLow=0.1, rhoHigh=0.5, verbose=False)

    lsparams = linesearchHZ.LinesearchHZ.defaultParameters()
    lsparams.sigma=0.1
    lsparams.delta=0.05
    lsparams.rename('linesearch')
    params.add(lsparams)

    linearparams = KrylovSolver.defaultParameters()
    linearparams.mu=1.
    linearparams.rename('linearsolver')
    params.add(linearparams)
    return params

  def __init__(self, params=None):
    self.params = self.defaultParameters()
    if not params is None: self.params.update(params)
    self.iteration_listeners = []
    self.linear_iteration_listeners = []
    self.x_update_listeners = []

  def addIterationListener(self,listener):
    """
    Add an object to be called after each iteration.
    """
    self.iteration_listeners.append(listener)

  def addLinearIterationListener(self,listener):
    """
    Add an object to be called after each iteration during solution
    of the linearized ill-posed problem.
    """
    self.linear_iteration_listeners.append(listener)

  def addXUpdateListener(self,listener):
    self.x_update_listeners.append(listener)

  ###################################################################################
  ##
  ## The methods up to the comment below are the ones that can or must be overridden by a subclass.
  ##
  ## Those implemented with a NotImpelmentedError must be overridden.

  def initialize(self,x,y,*args):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.

    The remaining arguments are passed directly from :func:`solve`, and might be used for determining the
    final stopping criterion.

    Returns vectors corresponding to the initial value of *x* and the desired value of *y=F(x)* along
    with the desired terminal discrepancy.
    """
    raise NotImplementedError()

  def discrepancy(self,x,y,r):
    """
    During the computation, a class-defined scalar called a discrepancy is maintained that 
    signifies in a class-defined way the difference between the currently computed value of F(x)
    and the desired value y.  Typically it will be
    
    .. math::         {\\rm discrepancy} = ||y-\calF(x)||_Y
        
    but there are scenarios where some other scalar is maintained.  During each linear step of
    the computations, a fraction of the discrepancy is eliminated.

    :param x: The current point in the domain.
    :param y: The desired value of F(x)
    :param r: The residual y-F(x)
    """
    return sqrt(abs(self.forwardProblem().rangeIP(r,r)))
    
  def forwardProblem(self):
    """
    Returns the :class:`NonlinearForwardProblem` that defines the maps F, T, T*, and the inner products on the 
    domain and range.
    """
    raise NotImplementedError()

  def temper_d(self,x,d,y,residual):
    pass

  def linearInverseSolve(self,x,y,r,discrepancyLin):
    """Internal method for solving the linearized inverse problems."""
    x0 = x.zero_like()
    
    forward_problem = self.forwardProblem()
    
    linSolver = BasicKrylovCGNE(forward_problem, params=self.params.linearsolver)
    
    for l in self.linear_iteration_listeners:
      linSolver.addIterationListener(l)
    
    (h,Th) = linSolver.solve(x0,r,discrepancyLin)
    return h
    
  def finalize(self,x,y):
    """
    Hook called at the end of :func:`solve`.  Gives the chance to massage the return values.
    """
    return (x,y)

  def iterationHook(self,count,x,Fx,y,d,r,Td):
    """
    Called during each iteration with the pertinent computations.  Handy for debugging and visualization.
    """
    for listener in self.iteration_listeners:
      listener(self,count,x,Fx,y,d,r,Td)

  def xUpdateHook(self,count,x,Fx,y,r):
    for listener in self.x_update_listeners:
      listener(self,count,x,Fx,y,r)

  ##
  ## End of methods to be overridden by a subclass.
  ##
  ########################################################################################################

  def solve( self, x0, y, *args ):
    """Main routine to solve the inverse problem F(x)=y.  Initial guess is x=x0.
    Extra arguments are passed to :func:`initialize`."""
    (x,y,targetDiscrepancy) = self.initialize(x0,y,*args)

    self.discrepancy_history=[]

    forward_problem = self.forwardProblem()
    params = self.params
    
    cg_reset = x.size()
    if( self.params.cg_reset != 0): cg_reset = self.params.cg_reset

    # Initial functional evalutation
    Fx = forward_problem.evalFandLinearize(x)

    # Prepare some storage
    Td = None

    residual = y.copy()
    residual.axpy(-1, Fx)

    discrepancy = self.discrepancy(x,y,residual);
    
    # The class that performs our linesearches.
    line_search = linesearchHZ.LinesearchHZ(params=self.params.linesearch)
    line_searchee = ForwardProblemLineSearchAdaptor(forward_problem)

    # Main loop
    count = 0
    theta = params.thetaMax;
    # try:
    for kkkkk in range(1):
      while True:
        self.discrepancy_history.append(discrepancy)

        if count > self.params.ITER_MAX:
          raise IterationCountFailure(self.params.ITER_MAX)
        count += 1

        if theta < self.params.thetaMin:
          raise NumericalFailure( 'Reached smallest trust region size.')

        # Error to correct:
        discrepancyLin = (1-theta)*discrepancy + theta*targetDiscrepancy
        msg('(%d) discrepancy: current %g linear goal:%g goal: %g\n---', count, discrepancy, discrepancyLin, targetDiscrepancy)

        if discrepancy <= self.params.mu*targetDiscrepancy:
          msg('done at iteration %d', count)
          break

        # try:
        d = self.linearInverseSolve(x,y,residual,discrepancyLin)
        # except Exception as e:
        #   theta *= self.params.kappaTrust
        #   msg('Exception during linear inverse solve:\n%s\nShriking theta to %g.',str(e),theta)
        #   continue

        # forward_problem.evalFandLinearize(x,out=Fx)
        # residual[:] = y
        # residual -= Fx
        Td = forward_problem.T(d,out=Td)
        Jp = -forward_problem.rangeIP(Td,residual)

        if Jp >= 0:
          theta *= self.params.kappaTrust
          msg('Model problem found an uphill direction.  Shrinking theta to %g.',theta)
          continue

          # % Sometimes in the initial stages, the linearization asks for an adjustment many orders of magnitude
          # % larger than the size of the coefficient gamma.  This is due to the very shallow derivatives
          # % in the coefficent function (and hence large derivatives in its inverse).  We've been 
          # % guaranteed that dh is pointing downhill, so scale it so that its size is on the order 
          # % of the size of the current gamma and try it out.  If this doesn't do a good job, we'll end
          # % up reducing theta later.
          # if(params.forceGammaPositive)

        self.temper_d(x,d,y,residual)

        self.iterationHook(count,x,Fx,y,d,residual,Td)

        # Do a linesearch in the determined direction.
        Phi = lambda t: line_searchee.eval(x,d,y,t)
        Jx = 0.5*forward_problem.rangeIP(residual,residual)
        line_search.search(Phi,Jx,Jp,1)
        if line_search.error():
          msg('Linesearch failed: %s. Shrinking theta.', line_search.errMsg );        
          theta *= self.params.kappaTrust
          continue

        discrepancyPrev = discrepancy
        t = line_search.value.t
        x.axpy(t,d)

        Fx.set(line_search.value.data)

        # forward_problem.evalFandLinearize(x,out=Fx,guess=Fx)
        residual.set(y)
        residual -= Fx
        discrepancy = self.discrepancy(x,y,residual);

        self.xUpdateHook(count,x,Fx,y,residual)
        
        # Check to see if we did a good job reducing the misfit (compared to the amount that we asked
        # the linearized problem to correct).
        rho = (discrepancy-discrepancyPrev)/(discrepancyLin-discrepancyPrev)
      
        # assert(rho>0)

        # Determine if the trust region requires any adjustment.
        if rho > self.params.rhoLow:
          # We have a good fit.
          if rho > self.params.rhoHigh:
            if theta < self.params.thetaMax:
              theta = min( theta/self.params.kappaTrust, self.params.thetaMax );
              msg( 'Very good decrease (rho=%g).  New theta %g', rho, theta );
            else:
              msg( 'Very good decrease (rho=%g).  Keeping theta %g', rho, theta );
          else:
            msg( 'Reasonable decrease (rho=%g). Keeping theta %g',rho, theta );
        else:
          # We have a lousy fit.  Try asking for less correction.
          theta = theta * params.kappaTrust;
          msg( 'Poor decrease (rho=%g) from %g to %g;', rho, discrepancyPrev, discrepancy  )
          msg( 'wanted %g.  New theta %g', discrepancyLin, theta );
    # except Exception as e:
    #   # Store the current x and y values in case they are interesting to the caller, then
    #   # re-raise the exception.
    #   self.finalState = self.finalize(x,Fx)
    #   raise e
      
    return self.finalize(x, Fx)

class BasicInvertIGN(InvertIGN):
  """Inversion of a forward problem using nonlinear conjugate gradient minimization
  and the Morozov discrepancy principle."""

  def __init__(self,forward_problem,params=None):
    InvertIGN.__init__(self,params=params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    """
    Returns the NonlinearForwardProblem that defines the inverse problem. 
    """
    return self.forward_problem

  def solve(self,x0,y,targetDiscrepancy):
    """
    Run the iterative method starting from the initial point *x0*.

    The third argument is the desired value of :math:`||y-T(x)||_Y`
    """
    return InvertIGN.solve(self,x0,y,targetDiscrepancy)

  def initialize(self,x0,y,targetDiscrepancy):
    """
    This method is a hook called at the beginning of a run.  It gives an opportunity for the class to 
    set up information needed to decide conditions for the final stopping criterion.

    It may be that the initial data 'x0' expresses the the initial data for the problem T(x)=y
    indirectly. Or it could be that x0 and y are expressed as dolfin.Function's rather than dolfind.GenericVectors.
    So initialize returns a triple of vectors (x0,y) which are possibly modified versions of the input data.

    The arguments \*args are passed directly from 'run'.
    """
    return (x0,y,targetDiscrepancy)


class InvertNLCG:
  """Implements solving an ill-posed minimization problem of the form
  
.. math::     J(x) = B(y-F(x),y-F(x))

where B is a symmetric positive semidefinite bilinear form and F is a (possibly nonlinear) function mapping between
two Hilbert spaces.

Minimization is done applying a nonlinear conjugate gradient algorithm until a stopping criterion is met.
"""


  ###################################################################################
  ##
  ## The methods up to the comment below are the ones that can or must be overridden by a subclass.
  ##
  ## Those implemented with a NotImpelmentedError must be overridden.

  def iterationHook(self,count,x,Fx,y,d,r,TStarR):
    """
    Called during each iteration with the pertinent computations.  Handy for debugging and visualization.
    """
    for listener in self.iteration_listeners:
      listener(self,count,x,Fx,y,d,r,TStarR)

  def xUpdateHook(self,count,x,Fx,y,r):
    for listener in self.x_update_listeners:
      listener(self,count,x,Fx,y,r)

  def stopConditionMet(self,count,x,Fx,y,r):
    """
    Determines if minimization should be halted (based, e.g. on a Morozov discrepancy principle)
    
    In: 
        * count: current iteration count
        * x:     point in domain of potential minimizer.
        * Fx:    value of nonlinear function at x
        * y:     desired value of F(x)
        * r:     current residual    
    """
    raise NotImplementedError()

  def initialize(self,x,y,*args):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.
    
    The remaining arguments are passed directly from 'solve', and might be used for determining the
    final stopping criterion.
    
    Returns dolfin vectors corresponding to the initial value of x and the desired value of y=F(x).    
    """
    raise NotImplementedError()

  def forwardProblem(self):
    """
    Returns the NonlinearForwardProblem that defines the maps F, T, T*, and the inner products on the 
    domain and range.
    """
    raise NotImplementedError()
    
  def finalize(self,x,y):
    """
    Hook called at the end of 'solve'.  Gives the chance to massage the return values.
    By default returns (x,y).
    """
    return (x,y)

  ## End of methods to be overridden by a subclass.
  ##
  ########################################################################################################

  @staticmethod
  def defaultParameters():
    """Parameters:
      
      * ITER_MAX: maximum iteration count
      * mu: scaling parameter for Morozov discrepency principle.
      * cg_reset: reset conjugate gradient method after this many iterations.  Zero for a default behaviour.
      * steepest_descent: Use steepest descent rather than full NLCG
      * linesearch: parameters to be passed to the linesearch algorithm
      * verbose: print out extra output
    """
    params = Parameters('InvNLCG', ITER_MAX=200, mu=1.1, cg_reset=0, deriv_eps=1, steepest_descent=False, verbose=False)
    lsparams = linesearchHZ.LinesearchHZ.defaultParameters()
    lsparams.sigma=0.1
    lsparams.delta=0.05
    lsparams.rename('linesearch')
    params.add(lsparams)
    return params
    
  def __init__(self, params=None):
    self.params = self.defaultParameters()
    if not params is None: self.params.update(params)
    self.iteration_listeners = []
    self.x_update_listeners = []

  def addIterationListener(self,listener):
    """
    Add an object to be called after each iteration.
    """
    self.iteration_listeners.append(listener)

  def addXUpdateListener(self,listener):
    """
    Add an object to be called after each new value of x is computed.
    """
    self.x_update_listeners.append(listener)

  def solve( self, x0, y, *args ):
    """Solve the inverse problem F(x)=y.  The initial estimate is x=x0.
    Any extra arguments are passed to :func:`initialize`.
    """
    (x,y) = self.initialize(x0,y,*args)

    # self.discrepancy_history=[]

    forward_problem = self.forwardProblem()
    cg_reset = x.size()
    if( self.params.cg_reset != 0): cg_reset = self.params.cg_reset

    # Initial functional evalutation
    if self.params.verbose: msg('initial evaluation')
    Fx = forward_problem.evalFandLinearize(x)

    residual = y.copy()
    residual.axpy(-1, Fx)

    Jx = 0.5*forward_problem.rangeIP(residual,residual)

    TStarR = forward_problem.TStar(residual)
    TStarRLast = TStarR.copy()

    # Compute our first guess for the descent direction.  
    d = TStarR.copy();
    Jp = -forward_problem.domainIP(TStarR,d)

    if self.params.verbose: msg('initial J %g Jp %g', Jx, Jp)

    # We need to give an initial guess for the stepsize in the linesearch routine
    # t0 = Jx/(1-Jp);
    t0 = Jx/(self.params.deriv_eps-Jp)

    # An analog of matlab's realmin
    realmin = np.finfo(np.double).tiny
    
    # The class that performs our linesearches.
    line_search = linesearchHZ.LinesearchHZ(params=self.params.linesearch)
    line_searchee = ForwardProblemLineSearchAdaptor(forward_problem)

    # Keep track of multiple line search failures.
    prevLineSearchFailed = False

    try:
      # Main loop
      count = 0
      while True:
        # self.discrepancy_history.append(sqrt(2*Jx))

        if count > self.params.ITER_MAX:
          raise IterationCountFailure(self.params.ITER_MAX)
        count += 1

        if self.stopConditionMet(count,x,Fx,y,residual):
          msg('done at iteration %d', count)
          break


        self.iterationHook(count,x,Fx,y,d,residual,TStarR)

        # Phi = lambda t: self.evalPhiAndPhiPrime(forward_problem,x,d,y,t)
        Phi = lambda t: line_searchee.eval(x,d,y,t)
        line_search.search(Phi,Jx,Jp,t0)
        if line_search.error():
          if prevLineSearchFailed:
            raise NumericalFailure( 'linesearch failed twice in a row: %s' % line_search.errMsg );
          else:
            msg('linesearch failed: %s, switching to steepest descent', line_search.errMsg );
            d = TStarR;
            t = 1/(self.params.deriv_eps-Jp);
            prevLineSearchFailed = True;
            continue      

        prevLineSearchFailed = False;

        t = line_search.value.t;
        x.axpy(t,d)
        TStarRLast.set(TStarR)

        Fx = forward_problem.evalFandLinearize(x,out=Fx,guess=Fx)
        residual.set(y)
        residual -= Fx

        self.xUpdateHook(count,x,Fx,y,residual)

        TStarR = forward_problem.TStar(residual,out=TStarR)

        Jx = 0.5*forward_problem.rangeIP(residual,residual)

        if self.params.steepest_descent:
          beta = 0
        else:
          # Polak-Ribiere
          beta = forward_problem.domainIP(TStarR,TStarR-TStarRLast)/forward_problem.domainIP(TStarRLast,TStarRLast)

          # If we have done more iterations than we have points, reset the conjugate gradient method.
          if count > cg_reset:
            beta = 0

        d *= beta
        d += TStarR

        JpLast = Jp
        Jp =  -forward_problem.domainIP(TStarR,d)
        if Jp >=0:
          if self.params.verbose:
            msg('found an uphill direction; resetting!');
          d.set(TStarR)
          Jp = -forward_problem.domainIP(TStarR,d);
          t0 = Jx/(self.params.deriv_eps-Jp);
        else:
          t0 =  t* min(10, JpLast/(Jp - realmin));
    except Exception as e:
      # Store the current x and y values in case they are interesting to the caller, then
      # re-raise the exception.
      import traceback
      traceback.print_exc()
      self.finalState = self.finalize(x,Fx)
      raise e

    return self.finalize(x, Fx)

class BasicInvertNLCG(InvertNLCG):
  """Inversion of a forward problem using nonlinear conjugate gradient minimization
  and the Morozov discrepancy principle."""
  
  def __init__(self,forward_problem,params=None):
    InvertNLCG.__init__(self,params=params)
    self.forward_problem = forward_problem
    
  def forwardProblem(self):
    """
    Returns the NonlinearForwardProblem that defines the inverse problem. 
    """
    return self.forward_problem

  def solve(self,x0,y,targetDiscrepancy):
    """
    Run the iterative method starting from the initial point *x0*.

    The third argument is the desired value of :math:`||y-T(x)||_Y`
    """
    return InvertNLCG.solve(self,x0,y,targetDiscrepancy)

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

  def stopConditionMet(self,count,x,Fx,y,r):
    """
    Determines if minimization should be halted (based, e.g. on a Morozov discrepancy principle)

    In: 
        * count: current iteration count
        * x:     point in domain of potential minimizer.
        * Fx:    value of nonlinear function at x
        * y:     desired value of F(x)
        * r:     current residual    
    """
    J = sqrt(abs(self.forward_problem.rangeIP(r,r)));
    Jgoal = self.params.mu*self.targetDiscrepancy
    
    msg('(%d) J=%g goal=%g',count,J,Jgoal)
    
    return J <= Jgoal

  
