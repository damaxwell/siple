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

# The following was copied from matplotlib, which copied a python recipe.
class Bunch:
    """
    Often we want to just collect a bunch of stuff together, naming each
    item of the bunch; a dictionary's OK for that, but a small do- nothing
    class is even handier, and prettier to use.  Whenever you want to
    group a few variables:

      >>> point = Bunch(datum=2, squared=4, coord=12)
      >>> point.datum
      By: Alex Martelli
      From: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52308
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def has_key(self,k):
      return self.__dict__.has_key(k)
      
    def update(self,**kwds):
      self.__dict__.update(**kwds)

    def __repr__(self):
        keys = self.__dict__.keys()
        return 'Bunch(%s)'%', '.join(['%s=%s'%(k,self.__dict__[k]) for k in keys])

class Parameters(Bunch):
  """This is a quick and dirty replacement for the dolfin.Parameters class
  that was previously used in 'sipl'."""
  def __init__(self, *args, **kwds):
      Bunch.__init__(self,**kwds)
      if len(args) > 0:
        self.name = args[0]

  def __repr__(self):
      keys = self.__dict__.keys()
      return 'Parameters(%s): [%s]'% (self.name,','.join(['%s=%s'%(k,self.__dict__[k]) for k in keys]) )

  def add(self,*args):
    if len(args)==0 or len(args)>2:
      'add expects a single Parameters argument or a key/value pair'
    if(len(args)==1):
      otherParams = args[0]
      if(isinstance(otherParams,Parameters)):
        self.__dict__[otherParams.name] = otherParams.copy()
      else:
        raise ValueError('add with one argument expects a Parameters, found %s' % otherParams)
    else:
      name = args[0]
      value = args[1]
      self.__dict__[name] = value
        
  def rename(self,name):
    self.name = name

  def update(self,other):
    if isinstance(other,dict):
      self.__dict__.update(other)
    elif isinstance(other,Parameters):
      self.__dict__.update(other.__dict__)
    else:
      raise ValueError()

  def copy(self):
    p=Parameters(self.name)
    p.__dict__.update(self.__dict__)
    return p