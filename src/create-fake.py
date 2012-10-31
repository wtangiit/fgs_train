#!/usr/bin/env python
'''This function creates an almost uninformative FGS training file
using only emission probabilities of .0001 for each of three hard-coded
stop codons and otherwise uninformative linear functions of GC content.'''

import sys, os
from optparse import OptionParser


if __name__ == '__main__':
  usage  = "usage: %prog -i <input sequence file> -o <output file>"
  parser = OptionParser(usage)
  parser.add_option("-i", "--input",  dest="input", default=None, help="Input sequence file.")
  parser.add_option("-o", "--output", dest="output", default=None, help="Output file.")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  parser.add_option("-r", "--reverse", dest="rev", action="store_true", default=False, help="Output rgene")
  
  parser.add_option("-n", "--noncoding", dest="non", action="store_true", default=False, help="Output noncoding")
  (opts, args) = parser.parse_args()
  rev = opts.rev
  non = opts.non
  for i in range(0, 45):
    gc = ( (float(i) +26) / 100) 
    print i+26
    if not non: 
     for j in range(0,6):
      for k in range(0,16):
       a = .5 -gc / 2
       t = .5 -gc / 2
       c =  gc / 2
       g =  gc / 2 
       if rev ==0: 
        if ( (j == 2 or j == 5) and k== 12 ):
          g=0.0001
          a=0.0001
          r = c+t
          c = c/(r +.0002)
          t = t/(r +.0002)
        if ( (j == 2 or j == 5) and k== 14 ):
          a=0.0001 
          r=c+g+t
          c = c/(r +.0001)
          g = g/(r +.0001)
          t = t/(r +.0001)
       if rev ==1:
        if ( (j == 2 or j == 5) and (k== 7 or k==13 or k==15) ):
          a=0.0001 
          r=c+g+t
          c = c/(r +.0001)
          g = g/(r +.0001)
          t = t/(r +.0001)
        if ( (j == 2 or j == 5) and k== 15 ):
          a=0.0001 
          r=c+g+t
          c = c/(r +.0001)
          g = g/(r +.0001)
          t = t/(r +.0001)
       print "%.04f\t%.04f\t%.04f\t%.04f"%( a, c, g, t )
    else:
       a = .5 -gc / 2
       t = .5 -gc / 2
       c =  gc / 2
       g =  gc / 2
       for i in range(0,4):
         print "%.04f\t%.04f\t%.04f\t%.04f"%( a, c, g, t )
#  if opts.verbose: sys.stdout.write("Done. \n")

