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

import inspect
import os.path
# This file contains functions helpful work working with python in an interactive session.

logFile = None
def logToFile(b):
  global logFile
  if isinstance(b,str):
    logFile = file(string,'w')
  else:
    if b:
      calling_file=inspect.getouterframes(inspect.currentframe())[1][1]
      log_file = os.path.splitext(calling_file)[0]+'.log'
      logFile=file(log_file,'w')
    else:
      if not logFile is None:
        logFile.close()
      logFile = None
      

def msg(s,*args):
  """
  Print a nicely formatted message.
  """
  global logFile
  framerec = inspect.stack()[1]
  calling_module = inspect.getmodule(framerec[0])
  f=inspect.getouterframes(inspect.currentframe())[1]
  caller_name = os.path.basename(f[1])
  if isinstance(f[3],str):
    caller_name += (':%s' % f[3])
  message = '%s: %s' % (caller_name, s % args)
  print( message )
  if not logFile is None:
    logFile.write(message)
    logFile.write('\n')
    logFile.flush()

def pause():
  """
  Halt computation until a key is pressed.
  """
  print ("Computation paused. Press any key to continue.")
  getch()
  print "Continued"

def beep():
  print '\a'


## The  following Getch fucntions are taken from {{{ http://code.activestate.com/recipes/134892/ (r2)
class _Getch:
    """Gets a single character from standard input.  Does not echo to the
screen."""
    def __init__(self):
        try:
            self.impl = _GetchWindows()
        except ImportError:
            self.impl = _GetchUnix()

    def __call__(self): return self.impl()

class _GetchUnix:
    def __init__(self):
        import sys, termios, fcntl, os

    def __call__(self):
        import sys, termios, fcntl, os
        fd = sys.stdin.fileno()
        oldattr = termios.tcgetattr(fd)
        oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
        c = 0
        try:        
          newattr = termios.tcgetattr(fd)
          newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
          termios.tcsetattr(fd, termios.TCSANOW, newattr)
          fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)
          
          while True:
            try:
              c=sys.stdin.read(1)
              break
            except IOError: pass

        finally:
          termios.tcsetattr(fd, termios.TCSADRAIN, oldattr)
          fcntl.fcntl(fd, fcntl.F_SETFL, oldflags)

        return c
      
class _GetchUnixOrig:
    def __init__(self):
        import tty, sys, termios

    def __call__(self):
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch


class _GetchWindows:
    def __init__(self):
        import msvcrt

    def __call__(self):
        import msvcrt
        return msvcrt.getch()


getch = _Getch()
## end of http://code.activestate.com/recipes/134892/ }}}
