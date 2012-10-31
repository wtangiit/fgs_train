#!/usr/bin/env python

'''do linear fit for generated gene.ct / rgene.ct file

take trained gene.ct or rgene.ct as input, generate a new trained file gene.fit or rgene.fit where
transition probability is the linear function of gc content. In the mean time, generate the 
plots showing linear regression.

./linear_fit.py  -i gene -w gene.ct
'''
import os
from optparse import OptionParser
import numpy as np
import scipy
from scipy.optimize import leastsq
from numpy import arange,array,ones
from pylab import plot,show
import matplotlib.pyplot as plt

digit2nt = {0: 'A', 1:'C', 2:'G', 3:'T'}

linestyle1 = ['ob', 'or', 'ok', 'og']
linestyle2 = ['-b', '-r', '-k', '-g',]

linfunc = lambda p,x:      p[0] + p[1]*x   
errfunc = lambda p,x,y,er : (y-linfunc(p,x) ) / er 

def number2dimer(number):
    digit1 = number / 4
    digit2 = number % 4
    nt1 = digit2nt[digit1]
    nt2 = digit2nt[digit2]
    return nt1+nt2    

def list_fit(counts, m, ij):
    ''' 
    counts              4x45 list containing counts for each of 4 nucleotides at each of 45 gc buckets
    m                   int  (state, M0, M1, etc)
    ij                  int  (dinucleotide AA, AC, etc)
       ''' 
    xi = arange(0,45)
    A = array([ xi, ones(45)])
    array_counts = np.array(counts, dtype=float)
    total_counts=array_counts.sum(axis=0)
    y = np.array(counts, dtype=float)
    ylist = array(y / total_counts, dtype=float)

    inv_weight = np.sqrt((ylist *(1-ylist) / total_counts))  # 4x45  same shape as ylist
    inv_weight = np.sqrt(( 1 / total_counts)) * np.ones([  ylist.shape[0] ,1]  )                # 4x45  same shape as ylist
    for i in  np.nonzero(total_counts==0): ylist[:,i]=0        # removes nan from divide
    for i in  np.nonzero(total_counts==0): inv_weight[:,i]=1 # assigns nonzero errorbars 
 #   print repr(ylist)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    from_dimer = number2dimer(ij)
    plt.title("M%d P(X|%s)" % (m, from_dimer))

    plots = []
    legend_labels = []
    
    parameter_list = []

    for k in range(4):
        y = ylist[k]
        p0 = [np.average(y), 0]
        print "p0 ", p0

        inverr=inv_weight[k]
        to_nt = digit2nt[k]

        w = scipy.optimize.leastsq(errfunc, p0, args=(xi,y,inverr),  maxfev=2000 )[0]
        parameter_list.append(w)
        
        line = w[0] + w[1]*xi # regression line
        plots.append(ax.plot(range(26, 71),line,linestyle2[k]))
        plots.append( [ax.errorbar(range(26, 71),y, yerr=inverr, fmt=linestyle1[k]  )[0] ])
        legend_labels.append(to_nt)
        
        print "plotting line M%s:P(%s|%s), slope=%.4f, intercept=%.4f" % (m, to_nt, from_dimer, w[1], w[0])
        
    if len(total_counts) > 0: #plot weight line
        ax2 = ax.twinx()
        plots.append(ax2.plot(range(26, 71), total_counts, '-xm', markeredgewidth=1))
        ax.set_ylim(0, 1)
#        ax2.set_ylim(0, 400000)
        ax2.set_ylabel('triplet count')
    
    
    ax.legend((plots[1][0], plots[3][0], plots[5][0], plots[7][0], ), ('A', 'C', 'G', 'T'), loc=0)

    ax.set_xlabel('GC content')
    ax.set_ylabel('probability')
    ax.grid(True)
    
    savefile = "M%d_%s.png" % (m, from_dimer)
    plt.savefig(savefile)
    
    return parameter_list
            
def parse_file(filename):
    infile = open(filename, "r")

    data_lists = [[[[] for k in range(4)] for ij in range(16)] for m in range(6)]# 6x16x4x[45]
        
    for g in range(45):
        line = infile.readline().strip('\n')
        for m in range(6):
            for ij in range(16):
                line = infile.readline().strip('\n')
                numbers = line.split()
                for k in range(4):
                    prob = numbers[k]
                    data_lists[m][ij][k].append(prob)
    return data_lists

def gen_fitted_file(linear_parameters):
    '''generated a new trained gene transition model with fitted linear functions'''
    
    fit_data_lists = [[[[] for ij in range(16)] for m in range(6)] for g in range(45)] #45x6x16x4
    
    for m in range(6):
        for ij in range(16):
            for k in range(4):
                w = linear_parameters[m][ij][k]
                slope = w[1]
                intercept = w[0]
                for g in range(45):
                    if g < 5:
                        xi = 5
                    else:
                        xi = g  
                    if g >40 :
                        xi = 40
                    else:
                        xi = g  

                    prob = slope * xi + intercept
                    if prob < 0.0001:
                        prob = 0.0001
                        
                    fit_data_lists[g][m][ij].append(prob)
    
    outfile = open("gene.fit", "w")
    
    for g in range(45):
        gc = g+26
        line = "%s\n" % gc
        outfile.write(line)
        for m in range(6):
            for ij in range(16):
                line = ""
                for k in range(4):
                    prob = fit_data_lists[g][m][ij][k]
                    line += "%.4f\t" % prob
                line.strip('\t')
                line += '\n'
                outfile.write(line)
    outfile.close()    
    
if __name__ == '__main__':
    usage  = "usage: %prog -w <training_gene.ct> -o <output gene>"
    parser = OptionParser(usage)
    parser.add_option("-w", "--weight",  dest="weight", type = "string", default=None, help="<input weighting data>")
    parser.add_option("-o", "--output",  dest="output", type = "string", default=None, help="<ouptut trained matrix>>")
    
    (opts, args) = parser.parse_args()
    
    weight_lists = parse_file(opts.weight)  #datalists array of 6x16x4x45  states * dinuc. * final nuc.  * gc bucket 
            
    out_handle = open(opts.output, "w")  

    linear_parameter = [[[]for ij in range(16)] for m in range(6)]
    
    for m in range(6):
        for ij in range(16):
            ylist = weight_lists[m][ij]    #ylist  4x45
            
            weight_line = []
            
            parameters = list_fit(ylist, m, ij )
            
            linear_parameter[m][ij] = parameters
                
    gen_fitted_file(linear_parameter)
    print "Done."
