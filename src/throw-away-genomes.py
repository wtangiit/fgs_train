#!/usr/bin/env python
# This script takes a list of genomes, gc contents, and genome sizes
# and outputs a subset of this list which is more homogeneous with
# respect to basepairs per gc content bucket.
# assumes GC content field is integers.

import sys, os, random
from optparse import OptionParser
import numpy as np

mingenomesize = 100000
genomesizeperbucket = 4000000

if __name__ == '__main__':
  usage  = "usage: %prog <input table> "
  parser = OptionParser(usage)
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  filename = args[0]
  if not (filename and os.path.isfile(filename) ):
    parser.error("Missing input file" )
  b = {}
  g = {}
  a1 = [0]* 100
  a2 = [0] * 100
  for line in open(filename):
     f=line.split()
     if f[0][0] != "#" :
      if int(f[2]) > mingenomesize:
       try :
          b[f[1]][f[0]] = f[2]
       except KeyError:
          b[f[1]] = {}
          b[f[1]][f[0]]  =  f[2] 
  genomebuckets = b.keys()
  for gc in sorted(b.keys()):
    totalsize = 0
    for genome in b[gc].keys():
      genomesize = int(b[gc][genome])
      totalsize += genomesize
      a1[int(gc)] = totalsize
    genomelist = b[gc].keys()
    random.shuffle(genomelist)
    sizesofar = 0
    for genome in genomelist:
       if sizesofar < genomesizeperbucket:
         print "select %d %s"%(int(gc), genome), 
         a2[int(gc)] += int(b[gc][genome])
       else:
         print "rejected %d %s"%(int(gc), genome) ,
       sizesofar += int(b[gc][genome])
       print sizesofar
  for i in range(0,100):
    print i, a1[i], a2[i]    
   
