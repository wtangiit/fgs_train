#!/usr/bin/env python
'''plot_hist.py
plot histogram plots to fit normal distribution for pwm scores.
input: list of pwm score values, one list per line 
output: the histogram fitting plots for each list of score values
 e.g. start.pwm etc, which is the output of calc_pwm.py with '-r' specified
'''

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import math
from optparse import OptionParser


def plot_hist(x, mu, sigma, filename, gc):
    # the histogram of the data
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    
    plt.xlabel('Scores')
    plt.ylabel('Probability')
    plt.title('%s gc=%s   mu=%.2f  sigma=%.2f'%(filename, gc, mu,sigma))
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    filename = "hist_%s_%s.png" % (filename, gc)  
    plt.savefig(filename)
    
    
    
if __name__ == '__main__':
    usage  = "usage: %prog -i <input file, e.g. start.pwm>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="input", type = "string", default=None, help="<input training data>")
    
    (opts, args) = parser.parse_args()
    
    infile = open(opts.input, "r")
    
    keyword = opts.input.split('.')[0]
    
    gc_content = 0
    for line in infile:
        line.strip('\n')
        numbers = line.split()
        if len(numbers) == 1:
            gc_content = int(numbers[0])
        else:
            numbers = [ float(x) for x in numbers ]
            mean = sum(numbers) / len(numbers)
            variance  = 0
            for number in numbers:
                variance += (number - mean)**2
            variance = variance / len(numbers)
            sigma = math.sqrt(variance)
            print "plotting %s %s numbers=%s mu=%s sigma=%s" % (keyword, gc_content, len(numbers), mean, sigma)
            plot_hist(numbers, mean, sigma, keyword, gc_content)
        