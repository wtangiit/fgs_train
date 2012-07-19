#!/usr/bin/env python

'''fgs_train.py: count gc content and bp of each genome in a folder

Usage:

count_gc_genomes.py -f input_folder

''' 


import os
from optparse import OptionParser
import glob

nt_list = ['A', 'C', 'G', 'T']

nt_dict = {'A':0, 'C': 1, 'G':2, 'T':3, 
           'a':0, 'c':1, 'g':2, 't':3}

digit2nt = {0: 'A', 1:'C', 2:'G', 3:'T'}

complement_dict = {'A':'T', 'C': 'G', 'G':'C', 'T':'A', 
           'a':'t', 'c':'g', 'g':'c', 't':'a'}

stop_codons = ["TAG", "TAA", "TGA"]
start1_codons = ["CTA", "TTA", "TCA"]
 
start_codons = ["ATG", "GTG", "TTG"]
stop1_codons = ["CAT", "CAC", "CAA"]

STRATIFY = False

MIN_GC_CONTENT = 26
MAX_GC_CONTENT = 70
NUM_STRATIFY = 45
NUM_M_STATE = 6
NUM_DIMER = 16
NUM_NT = 4


def trimer_to_int(triplet):
    '''return number by triplet'''
    t1 = nt_dict.get(triplet[0])
    t2 = nt_dict.get(triplet[1])
    t3 = nt_dict.get(triplet[2])
    
    if t1 >= 0 and t2 >=0 and t3 >=0:
        return t1 * 16 + t2 * 4 + t3
    else:
        return -1
    
def get_gc_content(sequence):
    '''return gc_content% of a given sequence'''
    gc_count = 0
    for ch in sequence:
        if ch in ['G', 'C', 'g', 'c']:
            gc_count += 1
#    print "%f\t%f\t%f"%(float(gc_count) / len(sequence), round(float(gc_count) / len(sequence), 2) * 100, int(round(float(gc_count) / len(sequence), 2) * 100)  ) 
    gc_content = int(float(gc_count) / len(sequence) * 100 + 0.5)
    #print "gc countent %s rounded to %s" % (float(gc_count) / len(sequence) * 100, gc_content)
    if gc_content < MIN_GC_CONTENT:
        gc_content = MIN_GC_CONTENT
    if gc_content > MAX_GC_CONTENT:
        gc_content = MAX_GC_CONTENT
    return gc_content

def parse_input_file(filename, noncoding=False):
    '''parse the input file'''
    infile = open(filename, "r")
    seq_lists = []
    
    ct_dict = {}
    
    total_count = 0
    for line in infile:
        line = line.strip('\n')
        line = line.strip('\r')
        splits = line.split('\t')
        if line[0] == '>':
            continue
        
        for ch in line:
            if ch in ct_dict.keys():
                ch = ch.upper()
                ct_dict[ch] += 1
            else:
                ct_dict[ch] = 1
            total_count += 1
            
    gc_ct = ct_dict.get('G', 0) + ct_dict.get('C', 0)
    gc_content = int(float(gc_ct) / total_count * 100 + 0.5)
    
    return gc_content, total_count
        
if __name__ == '__main__':
    usage  = "usage: %prog -i <input genome file> -f <input folder path>  (specify either -i or -f)"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="input", type = "string", default=None, help="Input genome file.")
    parser.add_option("-f", "--folder",  dest="folder", type = "string", default=None, help="folder containing input genome file.")
    
    (opts, args) = parser.parse_args()
        
    if opts.input:
        gc_content, total_ct =  parse_input_file(opts.input)
        print "%s gc=%s total_bp=%s" % (opts.input, gc_content, total_ct)
        
    if opts.folder:
        path = opts.folder + '/'
        for infile in glob.glob( os.path.join(path, '*.fna') ):
            gc_content, total_ct =  parse_input_file(infile)
            print "%s gc=%s total_bp=%s" % (infile, gc_content, total_ct)
        
    
    

        
    
