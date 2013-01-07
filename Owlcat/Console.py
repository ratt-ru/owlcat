# -*- coding: utf-8 -*-

import os
import sys
import time
#import curses

# This stuff doesn't work when running the script in the background!
# Dunno why, so disabling this for 1.2

## Get width of terminal
## I couldn't be bothered to get newline mode working properly when curses is active,
## so I used curses.wrapper() do init curses, get the width, then close curses.
#def get_width (scr):
  #global _width;
  #_width = scr.getmaxyx()[1] - 1;

## Since curses can fall over when invoked via ssh non-interactively (curse you curses!),
## protect this with an exception
#try:
  #curses.wrapper(get_width);
#except:
  #_width = 80;

try:
  _height,_width = map(int,os.popen('stty size', 'r').read().split())
except:
  _width = 132;

def timestamp (time_start,format="%H:%M:%S:"):
  return time.strftime(format,time.gmtime(time.time()-time_start));

class Reporter (object):
  """A Reporter is used to make progress reports on the console, with optional timestamps."""
  def __init__ (self,timestamp=False):
    self.time_start = timestamp and time.time();

  def overprint (self,message):
    return self.pprint(message+"\r");

  def pprint (self,message,newline=True):
    # print message
    if self.time_start:
      message = "%s %s"%(timestamp(self.time_start),message); 
    if message[-1] == "\r":
      endline = "\r";
    else:
      endline = "\n\r";
    # cut message at width of terminal, and pad with spaces
    sys.stdout.write("%-*.*s"%(_width,_width,message.strip()));
    sys.stdout.write(endline);
    sys.stdout.flush();

  def __call__ (self,*args):
    self.pprint(" ".join(args));
    
  
  