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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%
# %%  linesearchHZ
# %%
# %%  Finds an approximate minimizer of F(t) for t in [0,infinity) satisfying Wolfe conditions.
# %%
# %%  Algorithm from: Hager, W. and Zhang, H. CG DESCENT, a Conjugate Gradient Method with Guaranteed Descent
# %%                  Algorithm 851. ACM Transactions on Mathematical Software. 2006;32(1):113-137.
# %%
# %%  Usage:  status = linesearchHZ( F0, F0p, F, t, params )
# %%
# %%  In: 
# %%  F0      - F(0)
# %%  F0p     - F'(0)
# %%  F       - Function to minimize.  See requirements below.
# %%  t       - Initial guess for location of minimizer.
# %%  params  - Optional struct of control parameters.  See below.
# %%  
# %%  Out:
# %%  status  - Structure with the following fields:
# %%  code    - nonnegative integer with 0 indicating no error
# %%  msg     - if code>0, a descriptive error message about why the algorithm failed
# %%  val       - if code==0, structure containing information about the minimizer
# %%  val.t     - location of minimizer
# %%  val.F     - value of F(c.t)
# %%  val.Fp    - value of F'(c.t)
# %%  val.data  - additional data at minimizer.  See below.
# %%
# %%  The function F must have the following signature: function [f, fp, fdata] = F(t).
# %%  The f and fp are the values of the function and derivative at t.  For some functions,
# %%  there is expensive data that are computed along the way that might be needed by the end user.
# %%  The data field allows this excess data to be saved and returned.  In the end it will show up
# %%  in status.c.data.
# %%
# %%  The various parameters that control the algorithm are described in the above reference and
# %%  have the same name.  The principal ones you might want to change:
# %%
# %%  delta: Controls sufficent decrease Wolfe condition:  delta*F'(0) >= (F(t)-F(0))/t
# %%  sigma: Controls sufficent shallowness Wolfe condition:  F'(t) >= sigma*F'(0)
# %%  rho:   Expansion factor for initial bracket search (bracket expands by multiples of rho)
# %%  nsecant: Maximum number of outermost loops
# %%  nshrink: Maximum number of interval shrinks or expansions in a single loop
# %%  verbose: Print lots of messages to track the algorithm
# %%  debug:   Do extra computations to verify the state is consistant as the algorithm progresses.
# %%
# %%  
# 

from siple.params import Bunch, Parameters
from siple.reporting import msg, pause
import numpy

class LinesearchHZ:
  
  @staticmethod
  def defaultParameters():
    return Parameters( 'linesearchHZ', delta=.1, sigma=.9, epsilon=0, theta=.5, gamma=.66,  rho=5,
                              nsecant=50, nshrink=50, verbose=False, debug=True );

  def __init__(self,params=None):
    self.params = self.defaultParameters()
    if not (params is None): self.params.update(params)

  def error(self):
    return self.code > 0
    
  def ezsearch(self,F,t0=None):
    self.F = F
    z = self.eval(0)
    if t0 is None:
      t0 = 1./(1.-z.F0p);
    return self.search(F,z.F,z.Fp,t0)
    
  def search(self,F,F0,F0p,t0):
    self.code = -1
    self.errMsg = 'no error'
    self.F = F

    params = self.params

    z = Bunch(F=F0,Fp=F0p,t=0,data=None)
    assert F0p <= 0

    # % Set up constants for checking Wolfe conditions.
    self.wolfe_lo = params.sigma*z.Fp;
    self.wolfe_hi = params.delta*z.Fp;
    self.awolfe_hi = (2*params.delta-1)*z.Fp;
    self.fpert = z.F + params.epsilon;
    self.f0 = z.F;

    if params.verbose: msg('starting at z=%g (%g,%g)', z.t, z.F, z.Fp )

    while True:
      c = self.eval(t0)
      if not numpy.isnan(c.F):
        break
      msg('Hit a NaN in initial evaluation at t=%g',t0)        
      t0 *= 0.5
    
    if params.verbose: msg('initial guess c=%g (%g,%g)', c.t, c.F, c.Fp )

    if self.wolfe(c):
      if params.verbose: msg('done at init')
      self.setDone(c)
      return

    (aj,bj) = self.bracket(z,c)
    if params.verbose: msg('initial bracket %g %g',aj.t,bj.t)
   
    if self.code >= 0:
      self.doneMsg('initial bracket')
      return

    if params.debug: self.verifyBracket(aj,bj)

    count = 0;

    while True:
      count += 1;

      if count> params.nsecant:
        self.setError('too many bisections in main loop')
        return

      (a,b) = self.secantsq(aj,bj);
      if params.verbose: msg('secantsq a %g b %g', a.t, b.t)
      if params.verbose: self.printBracket(a,b)
      if self.code >= 0:
        self.doneMsg('secant');
        return
      
      if  (b.t-a.t) > params.gamma*(bj.t-aj.t):
        (a,b) = self.update(a, b, (a.t+b.t)/2 );      
        if params.verbose: msg('update to a %g b %g', aj.t, bj.t)
        if params.verbose: self.printBracket(a,b)
        if self.code >= 0:
          self.doneMsg( 'bisect');
          return
      aj = a
      bj = b

  def printBracket(self,a,b):
    msg('a %g b %g f(a) %g fp(a) %g f(b) %g fp(b) %g fpert %g', a.t, b.t, a.F, a.Fp, b.F, b.Fp, self.fpert )

  def doneMsg(self,where):
    if self.code > 0:
      msg('done at %s with error status: %s', where, self.errMsg);
    else:
      if self.params.verbose: msg('done at %s with val=%g (%g, %g)', where, self.value.t, self.value.F, self.value.Fp );

  def verifyBracket(self,a,b):
    good = (a.Fp<=0) and (b.Fp >=0) and (a.F<= self.fpert);
    if not good:
      msg( 'bracket inconsistant: a %g b %g f(a) %g fp(a) %g f(b) %g fp(b) %g fpert %g', a.t, b.t, a.F, a.Fp, b.F, b.Fp, self.fpert )
      pause()
    if(a.t>=b.t):
      msg('bracket not a bracket (a>=b): a %g b %g f(a) %g fp(a) %g f(b) %g fp(b) %g fpert %g', a.t, b.t, a.F, a.Fp, b.F, b.Fp, self.fpert )

  def setDone(self, c):
    self.code = 0
    self.value = c
    
  def setError(self,msg):
    self.code = 1
    self.errMsg = msg

  def update(self, a, b, ct ):
    abar = a
    bbar = b
    
    params = self.params
    
    if params.verbose: msg('update %g %g %g', a.t, b.t, ct);
    if  (ct<=a.t) or (ct>=b.t):
      if params.verbose: msg('midpoint out of interval')
      return (abar,bbar)

    c = self.eval(ct)

    if self.wolfe(c):
      self.setDone(c)
      return (abar,bbar)

    if c.Fp >= 0:
      if params.verbose: msg('midpoint with non-negative slope. Becomes b.')
      abar = a;
      bbar = c;
      if params.debug: self.verifyBracket(abar,bbar)
      return (abar,bbar)

    if c.F <= self.fpert:
      if params.verbose: msg('midpoint with negative slope, small value. Becomes a.')
      abar = c;
      bbar = b;
      if params.debug: self.verifyBracket(abar,bbar)
      return (abar,bbar)

    if params.verbose: msg('midpoint with negative slope, large value. Shrinking to left.')
    (abar,bbar) = self.ushrink( a, c );
    if params.debug: self.verifyBracket(abar,bbar)
    
    return (abar,bbar)

  def ushrink(self,a,b):
    abar = a;
    bbar = b;

    count = 0;
    while True:
      count += 1;
      
      if self.params.verbose:
        msg('in ushrink')
        self.printBracket(abar,bbar)
      if count > self.params.nshrink:
        self.setError('too many contractions in ushrink')
        return (abar,bbar)
    
      d=self.eval( (1-self.params.theta)*abar.t+self.params.theta*bbar.t );
      if self.wolfe(d):
        self.setDone(d)
        return (abar,bbar)

      if d.Fp>=0:
        bbar = d;
        return (abar,bbar)
    
      if d.F <= self.fpert:
        abar=d;
      else:
        bbar=d;

  def plotInterval(self,a,b,N=20):
    from matplotlib import pyplot as pp
    import numpy as np
    T=np.linspace(a.t,b.t,N)
    FT=[]
    FpT=[]
    for t in T:
      c=self.eval(t)
      FT.append(c.F)
      FpT.append(c.Fp)
    pp.subplot(1,2,1)
    pp.plot(T,np.array(FT))
    pp.subplot(1,2,2)
    pp.plot(T,np.array(FpT))
    pp.draw()

  def secant(self,a,b):
    # % What if a'=b'?  We'll generate a +/-Inf, which will subsequently test as being out
    # % of any interval when 'update' is subsequently called.  So this seems safe.
    
    if self.params.verbose: msg( 'secant: a %g fp(a) %4.8g b %g fp(b) %4.8g',a.t,a.Fp, b.t, b.Fp)
    if(a.t==b.t):
      msg('a=b, inconcievable!')
    if -a.Fp <= b.Fp:
      return a.t-(a.t-b.t)*(a.Fp/(a.Fp-b.Fp));
    else:
      return b.t-(a.t-b.t)*((b.Fp)/(a.Fp-b.Fp));

  def secantsq(self,a,b):
    ct = self.secant(a,b)
    if self.params.verbose: msg('first secant to %g', ct)
    (A,B) = self.update(a,b,ct)
    if self.code >= 0:
      return (A,B)

    if B.t == ct:
      ct2 = self.secant(b,B);
      if self.params.verbose: msg('second secant on left half A %g B %g with c=%g',A.t, B.t, ct2)
      (abar,bbar) = self.update(A,B,ct2)
    elif A.t == ct:
      ct2 = self.secant(a,A);
      if self.params.verbose: msg('second secant on right half A %g B %g with c=%g',A.t, B.t, ct2)
      (abar,bbar) = self.update(A,B,ct2)
    else:
      if self.params.verbose: msg('first secant gave a shrink in update. Keeping A %g B %g',A.t, B.t)
      abar = A; bbar = B
    
    return (abar,bbar)


  def bracket( self, z, c ):
    a = z
    b = c
    
    count = 0
    while True:
      if count > self.params.nshrink:
        self.setError('Too many expansions in bracket')
        return (a,b)
      count += 1

      if b.Fp >= 0:
        if self.params.verbose: msg('initial bracket ends with expansion: b has positive slope')
        return (a,b)
    
      if b.F > self.fpert:
        if self.params.verbose: msg('initial bracket contraction')
        return self.ushrink(a,b);

      if self.params.verbose: msg('initial bracket expanding')
      a = b;
      rho = self.params.rho
      while True:
        if count > self.params.nshrink:
          self.setError('Unable to find a valid input')
          return (a,b)
        c = self.eval(rho*b.t)
        if not numpy.isnan(c.F):
          b = c
          break
        msg('Hit a NaN at t=%g',rho*b.t)
        rho*=0.5
        count += 1

      if self.wolfe(b):
        #msg('decrease %g slope %g f0 %g fpert %g', b.F-params.f0, b.t*params.wolfe_hi, params.f0, params.fpert)
        self.setDone(b);
        return(a,b)

  def wolfe(self,c):
    if self.params.verbose: msg('checking wolfe of c=%g (%g,%g)',c.t,c.F,c.Fp)

    if c.Fp >= self.wolfe_lo:
      if (c.F-self.f0) <= c.t*self.wolfe_hi:
        return True

      if self.params.verbose: msg('failed sufficient decrease' )

      # % if (  (c.F <= params.fpert) && (c.Fp <= params.awolfe_hi ) )
      # %   msg('met awolfe')
      # %   met = true;
      # %   return;
      # % end
      # if params.verbose: msg('failed awolfe sufficient decrease' )
    else:
      if self.params.verbose: msg('failed slope flatness')

    return False

  def eval(self,t):
    c = Bunch(F=0,Fp=0,data=None,t=t)
    (c.F,c.Fp,c.data) = self.F(t)
    c.F = float(c.F)
    c.Fp = float(c.Fp)
    return c

if __name__ == '__main__':

  lsParams = Parameters('tmp', verbose=True, debug=True)
  ls = LinesearchHZ(params=lsParams)
  F = lambda t: (-t*(1-t), -1+2*t,None)
  ls.ezsearch(F,5)
  if ls.error():
    print ls.errMsg
  else:
    v = ls.value
    print 'minimum of %g at t=%g' % (v.F,v.t)
