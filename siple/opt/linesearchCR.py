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

from numpy import isinf, isnan, sqrt, any, isreal, real, nan, inf
from siple.reporting import msg
from siple.params import Bunch, Parameters

#This program is distributed WITHOUT ANY WARRANTY; without even the implied 
#warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
#  This is a modified version of:
#
#
#This file contains a Python version of Carl Rasmussen's Matlab-function 
#minimize.m
#
#minimize.m is copyright (C) 1999 - 2006, Carl Edward Rasmussen.
#Python adaptation by Roland Memisevic 2008.
#
#
#The following is the original copyright notice that comes with the 
#function minimize.m
#(from http://www.kyb.tuebingen.mpg.de/bs/people/carl/code/minimize/Copyright):
#
#
#"(C) Copyright 1999 - 2006, Carl Edward Rasmussen
#
#Permission is granted for anyone to copy, use, or modify these
#programs and accompanying documents for purposes of research or
#education, provided this copyright notice is retained, and note is
#made of any changes that have been made.
#
#These programs and documents are distributed without any warranty,
#express or implied.  As the programs were written for research
#purposes only, they have not been tested to the degree that would be
#advisable in any important application.  All use of these programs is
#entirely at the user's own risk."

class LinesearchCR:
  
  @staticmethod
  def defaultParameters():
    return Parameters( 'linesearchCR',
    INT= 0.1, # don't reevaluate within 0.1 of the limit of the current bracket
    EXT = 3.0,              # extrapolate maximum 3 times the current step-size
    MAX = 20,                     # max 20 function evaluations per line search
    RATIO = 10,                                   # maximum allowed slope ratio
    SIG = 0.1,
    verbose=False ) # RHO = SIG/2 ) # SIG and RHO are the constants controlling the Wolfe-

  def __init__(self,params=None):
    self.params = self.defaultParameters()
    if not (params is None): self.params.update(params)

  def error(self):
    return self.code > 0

  def search(self,f,f0,df0,t0=None):
    INT = self.params.INT; # don't reevaluate within 0.1 of the limit of the current bracket
    EXT = self.params.EXT;              # extrapolate maximum 3 times the current step-size
    MAX = self.params.MAX;                     # max 20 function evaluations per line search
    RATIO = self.params.RATIO;                                   # maximum allowed slope ratio
    SIG = self.params.SIG;
    RHO = SIG/2; 
    SMALL = 10.**-16                    #minimize.m uses matlab's realmin 

    # SIG and RHO are the constants controlling the Wolfe-
    #Powell conditions. SIG is the maximum allowed absolute ratio between
    #previous and new slopes (derivatives in the search direction), thus setting
    #SIG to low (positive) values forces higher precision in the line-searches.
    #RHO is the minimum allowed fraction of the expected (from the slope at the
    #initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
    #Tuning of SIG (depending on the nature of the function to be optimized) may
    #speed up the minimization; it is probably not worth playing much with RHO.

    d0 = df0;
    fdata0 = None
    if t0 is None:
      t0 = 1/(1.-d0)
    x3 = t0                                   # initial step is red/(|s|+1)
    X0 = 0; F0 = f0; dF0 = df0; Fdata0=fdata0              # make a copy of current values
    M = MAX

    x2 = 0; f2 = f0; d2 = d0; 
    while True:                      # keep extrapolating as long as necessary
      # x2 = 0; f2 = f0; d2 = d0; 
      f3 = f0; df3 = df0; fdata3 = fdata0;
      success = 0
      while (not success) and (M > 0):
        try:
          M = M - 1
          (f3, df3, fdata3) = f(x3)
          if isnan(f3) or isinf(f3) or any(isnan(df3)+isinf(df3)):
            raise Exception('nan')
          success = 1
        except:                    # catch any error which occured in f
          if self.params.verbose: msg('error on extrapolate. shrinking %g to %g', x3, (x2+x3)/2)
          x3 = (x2+x3)/2                       # bisect and try again

      if f3 < F0:
        X0 = x3; F0 = f3; dF0 = df3; Fdata0=fdata3   # keep best values
      d3 = df3                                         # new slope
      if d3 > SIG*d0 or f3 > f0+x3*RHO*d0 or M == 0:  break # are we done extrapolati

      x1 = x2; f1 = f2; d1 = d2                 # move point 2 to point 1
      x2 = x3; f2 = f3; d2 = d3                 # move point 3 to point 2
      A = 6*(f1-f2)+3*(d2+d1)*(x2-x1)          # make cubic extrapolation
      B = 3*(f2-f1)-(2*d1+d2)*(x2-x1)
      Z = B+sqrt(complex(B*B-A*d1*(x2-x1)))
      if Z != 0.0:
          x3 = x1-d1*(x2-x1)**2/Z              # num. error possible, ok!
      else: 
          x3 = inf
      if (not isreal(x3)) or isnan(x3) or isinf(x3) or (x3 < 0): 
                                                 # num prob | wrong sign?
          x3 = x2*EXT                        # extrapolate maximum amount
      elif x3 > x2*EXT:           # new point beyond extrapolation limit?
          x3 = x2*EXT                        # extrapolate maximum amount
      elif x3 < x2+INT*(x2-x1):  # new point too close to previous point?
          x3 = x2+INT*(x2-x1)
      x3 = real(x3)
      msg('extrapolating: x1 %g d1 %g x2 %g d2 %g x3 %g',x1,d1,x2,d3,x3)


    while (abs(d3) > -SIG*d0 or f3 > f0+x3*RHO*d0) and M > 0:  
      if (d3 > 0) or (f3 > f0+x3*RHO*d0):            # choose subinterval
        x4 = x3; f4 = f3; d4 = d3             # move point 3 to point 4
      else:
        x2 = x3; f2 = f3; d2 = d3             # move point 3 to point 2
      if self.params.verbose: msg('interpolating x2 %g x4 %g f2 %g f4 %g wolfef %g d2 %g d4 %g wolfed %g ',x2,x4,f2,f4,f0+x3*RHO*d0,d2,d4,-SIG*d0)

      if f4 > f0:           
        x3 = x2-(0.5*d2*(x4-x2)**2)/(f4-f2-d2*(x4-x2)) # quadratic interpolation
      else:
        A = 6*(f2-f4)/(x4-x2)+3*(d4+d2)           # cubic interpolation
        B = 3*(f4-f2)-(2*d2+d4)*(x4-x2)
        if A != 0:
          x3=x2+(sqrt(B*B-A*d2*(x4-x2)**2)-B)/A # num. error possible, ok!
        else:
          x3 = inf
      if isnan(x3) or isinf(x3):
        x3 = (x2+x4)/2      # if we had a numerical problem then bisect
      x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2))  # don't accept too close
      (f3, df3, fdata3) = f(x3);
      d3 =df3;
      M = M - 1;                      # count epochs
      if f3 < F0:
        X0 = x3; F0 = f3; dF0 = df3; Fdata0 = fdata3              # keep best values

    if (abs(d3) < -SIG*d0) and (f3 < f0+x3*RHO*d0):          # if line search succeeded
      self.code = 0
      self.value = Bunch(F=f3,Fp=d3,t=x3,data=fdata3)
      self.errMsg = ""
    else:
      self.code = 1
      if M == 0:
        self.errMsg = 'Too many function evaluations (>%d)' % MAX
      else:
        self.errMsg = 'unknown error';
      self.value = Bunch(F=f0,Fp=dF0,t=X0,data=Fdata0)
