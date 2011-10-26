############################################################################
#
#  This file is a part of siple.
#
#  Copyright 2010-2011 David Maxwell
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
  """Indicates if |siple| output should be saved to a log file.  If *b* is a file name, output
  is saved there.  Otherwise *b* should be a boolean.  If it is :data:`True`, the log file name
  is taken from the name of the calling python file (with '.py' replaced with '.log').
  If *b* is :data:`False` then file logging is turned off.
  """
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

kPRATTLE=4
kDEBUG=3
kMESSAGE=2
kWARNING=1
kERROR=0
kSeverityDesc = [ "Error", "Warning", "Message", "DEBUG", "BLAHBLAH"]

# Default loggers
def logprint(message,severity):
  print message
  if severity <= kWARNING:
    beep()

def logfile(message,severity):
  global logFile
  if not logFile is None:
    logFile.write(message)
    logFile.write('\n')
    logFile.flush()

loggers = [ logprint, logfile ]

def clear_loggers():
  global loggers
  loggers = []

def add_logger(logger):
  global loggers
  loggers.append(logger)

def _format_message(caller_level,severity,s,*args):
  framerec = inspect.stack()[1+caller_level]
  calling_module = inspect.getmodule(framerec[0])
  f=inspect.getouterframes(inspect.currentframe())[1+caller_level]
  caller_name = os.path.basename(f[1])
  if isinstance(f[3],str):
    caller_name += (':%s' % f[3])
  if severity == kMESSAGE:
    sevdesc = ""
  else:
    sevdesc = "(%s)" % kSeverityDesc[severity]
  
  message = '%s:%s %s' % (caller_name, sevdesc, s % args)
  return message

format_message = _format_message
def set_message_formatter(formatter):
  global format_message
  format_message = formatter

def msg(s,*args):
  """
  Print a nicely formatted message.
  """
  global loggers
  caller_level = 1
  severity=kMESSAGE
  message = format_message(caller_level,severity,s,*args)
  for l in loggers:
    l(message,severity)

def prattle(s,*args):
  """
  Print a nicely formatted message.
  """
  global loggers
  caller_level=1
  severity=kPRATTLE
  message = _format_message(caller_level,severity,*args)
  for l in loggers:
    l(message,severity)

def debug(s,*args):
  global loggers
  caller_level=1
  severity=kDEBUG
  message = _format_message(caller_level,severity,*args)
  for l in loggers:
    l(message,severity)

def pause(message_in="Computation paused. Press any key to continue.",
          message_out='Continued'):
  pause_callback(message_in,message_out)

def set_pause_callback(callback):
  global pause_callback
  pause_callback = callback

def std_pause(message_in=None,message_out=None):
  """
  Halt computation until a key is pressed.
  """
  if not message_in is None:
    print(message_in)
  getch()
  if not message_out is None:
    print message_out

pause_callback = std_pause



def in_ipython():
    "Determines if the python interpreter is ipython."
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

def endpause(message_in='Script done. Press any key to continue',message_out=None,
             python_only=True):
  """
  Wait for a keypress after displaying a helpful message. Typically used at 
  the end of a script.  If python_only=True, the message and wait are skipped
  when running under ipython.
  """
  if (not python_only) or (not in_ipython()):
    pause(message_in,message_out)

def beep():
  "Beeps."
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
