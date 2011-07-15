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

class NumericalFailure(Exception):
  """
  Class for exceptoions indicating a failure of a numerical method.
  """
  pass
  
class IterationCountFailure(NumericalFailure):
  """
  A standard error: an interation count was exceeded.
  """
  def __init__(self,ITER_MAX):
    NumericalFailure.__init__(self,'Iteration count %d exceeded' % ITER_MAX)
    