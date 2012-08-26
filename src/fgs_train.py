#!/usr/bin/env python

'''fgs_train.py: produce needed training data set for running FragGeneScan HMM model
file can be trained: gene, rgene, noncoding, start, stop, start1, stop1.

Usage:

  fgs_train./py [-i <input sequence file>  | -d <input directory>] [-n <input noncoding file>] [-g]"

either a sequence file or a directory should be specified (can not specify both)
if specify an input sequence file, the single input file is used in training
if specify a directory containing sequence files, all the files within that directory are used in training

it will stratify genes by gc content (from 26% to 70%) unless setting '-g False' or '--gc_content=False' 
the output file also contains 45 groups of probability data, only difference 
is that every group has the same data.
''' 


import os
from optparse import OptionParser

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
    linenumber = 0    
    for line in infile:
        line = line.strip('\n')
        line = line.strip('\r')
        splits = line.split('\t')
        linenumber += 1
        if len(splits) == 3:
            seq_len = len(splits[2])
            if seq_len > 0:
                if False and not noncoding and seq_len < 123:
                    print "**********************************************"
                    print "Warning!!:the input data contains invalid data, line %d, sequence length=%s" % (linenumber, seq_len)
                    print "The invalid sequence is thrown out to continue, but replacing input data and re-training is suggested." 
                    print "**********************************************"
                    continue
                seq_lists.append(splits[2])
    return seq_lists

def parse_input_dir(input_dir):
    '''parse the input files in the given input_dir'''
    listing = os.listdir(input_dir)
    gene_list = []
    noncoding_list = [] 
    for infile in listing:
        if infile.find("fwd") > 0:
            seq_list = parse_input_file(infile)
            gene_list.extend(seq_list)
        if infile.find("noncoding") > 0:
            seq_list = parse_input_file(infile, noncoding=False)
            noncoding_list.extend(seq_list)
    return gene_list, noncoding_list
    
def train_gene_transition(seq_list, output_file):
    '''train transition probability of matching states'''
    
    e_M_counts = [[[[0 for i in range(4)] for j in range(16)] for m in range(6)] for g in range(NUM_STRATIFY) ]
    
    if output_file == "gene":
        stop_codon_list = stop_codons
    if output_file == "rgene":
        stop_codon_list = start1_codons
        
    for seq in seq_list:
        
        if STRATIFY:
            gc_content = get_gc_content(seq)
            
        else:
            gc_content = MIN_GC_CONTENT        
        
#        print "gc_content=", gc_content
        
        for i in range(60, len(seq)-63):  #iterate coding NTs
            m = i % 6  #m = 0..5 representing M state from 1..6
            to = nt_dict.get(seq[i], -1)
            from0 = nt_dict.get(seq[i-2], -1)
            from1 = nt_dict.get(seq[i-1], -1)
            if from0 >= 0 and from1 >= 0:
                from2 = from0 * 4 + from1
                e_M_counts[gc_content-MIN_GC_CONTENT][m][from2][to] += 1
            
    gene_file = open(output_file, "w")
    ct_file = open(output_file+".ct", "w")
    
    for gc in range(MIN_GC_CONTENT, MAX_GC_CONTENT + 1):
        line = "%s\n" % gc
        gene_file.write(line)
        ct_file.write(line)
        
        if STRATIFY:
            k = gc
        else:
            k = MIN_GC_CONTENT  
        
        for m in range(6):
        #print "position=", m+1
            for j in range(16):
                total_ct = sum(e_M_counts[k - MIN_GC_CONTENT][m][j])
                #  print dimer_list[j],
                line = ""
                line_ct = ""
                for i in range(4):
                    if total_ct > 0:
                        codon = "%s%s%s" % (digit2nt[j/4], digit2nt[j%4], digit2nt[i])
                        
                        ct = e_M_counts[k-MIN_GC_CONTENT][m][j][i]
                        prob = round(float(ct) / total_ct, 4)
                        
                        if prob < 0.001:
                            prob = 0.001
                            if codon in stop_codon_list and m in [2,5]:
                                prob = 0.0001
                                
                        if codon in stop_codon_list and m in [2,5]:
                            if ct > 0:
                                print "inframe stop (start1) codon found %s: gc=%s, m=%s, count=%s/%s, prob=%f" % (codon, gc, m, ct, total_ct, float(ct)/total_ct) 
             
                    else:
                        prob = 0.0001
                    line += str(prob)
                    line += '\t'
                    
                    line_ct += str(total_ct)
                    line_ct += '\t'                    
                    
                line = line.strip('\t')
                line += '\n'
                gene_file.write(line)
                
                line_ct = line_ct.strip('\t')
                line_ct += '\n'
                ct_file.write(line_ct)
    
    gene_file.close()
    ct_file.close()
    print "output file produced: %s" % output_file 
     
def train_gene_transition_two_way(seq_list):
    '''train gene transition probability files for two ways'''
    train_gene_transition(seq_list, "gene")
    rc_req_list = []
    for seq in seq_list:
        rc_req_list.append(get_reverse_complement(seq))
    train_gene_transition(rc_req_list, "rgene") 
    
def get_start_stop_subseq(seq, key):
    '''return the subsequence to check for training start/stop/start1/stop1'''
    
    rc_seq = get_reverse_complement(seq)
    
    if key=="start":
        subseq = seq[30:93]
    elif key=="stop":
        subseq = seq[-123:-60]
    elif key=="start1":
        subseq = rc_seq[-93:-30]
    elif key=="stop1":
        subseq = rc_seq[60:123]
    else:
        subseq = seq
    return subseq       
    
    
def train_start_stop_adjacent_prob(seq_list):
    '''train start, stop, start1, stop1, stratify by gene GC content'''
    
    prob_counts_dict = {"start": [[[0 for i in range(64)] for j in range(61)] for g in range(NUM_STRATIFY)],
                        "stop" : [[[0 for i in range(64)] for j in range(61)] for g in range(NUM_STRATIFY)],
                        "start1": [[[0 for i in range(64)] for j in range(61)] for g in range(NUM_STRATIFY)],
                        "stop1": [[[0 for i in range(64)] for j in range(61)] for g in range(NUM_STRATIFY)]}
                        
    for seq in seq_list:
        
        if STRATIFY:
            gc_content = get_gc_content(seq)
        else:
            gc_content = MIN_GC_CONTENT        
            #os.stderr.write("%d\n"%gc_content); 
        for key in prob_counts_dict.keys():
            subseq = get_start_stop_subseq(seq, key)
            for i in range(61):
                s_triplet = subseq[i:i+3]
                index = trimer_to_int(s_triplet)
                prob_counts_dict[key][gc_content-MIN_GC_CONTENT][i][index] += 1
    
    for key in prob_counts_dict.keys():
        write_start_stop_file(key, prob_counts_dict[key])
              
def write_start_stop_file(filename, prob_counts):
    '''write start stop prob into output files, with gc'''
    outfile = open(filename, "w")
    for gc in range(MIN_GC_CONTENT, MAX_GC_CONTENT + 1):
        line = "%s\n" % gc
        outfile.write(line)
        
        if STRATIFY:
            k = gc
        else:
            k = MIN_GC_CONTENT
        
        for i in range(61):
            line = "";
            total_ct = sum(prob_counts[k - MIN_GC_CONTENT][i])
            for j in range(64):
                if total_ct > 0:
                    prob = round(float(prob_counts[k - MIN_GC_CONTENT][i][j] + 1) / (total_ct + 1), 6)
                else:
                    prob = 0.000001
                if prob < 0.000001:
                    prob = 0.000001
                line += str(prob)
                line += '\t'
            line = line.strip('\t')
            line += ('\n')
            outfile.write(line)
    outfile.close()
    print "output file produced: %s" % filename
            
def get_reverse_complement(seq):
    '''return the reverse complement of the given sequence'''
    seq = seq[::-1]
    rseq= ""
    for ch in seq:
        try:
            rseq += complement_dict[ch]
        except KeyError:
            rseq += "N"
    return rseq            

def train_non_coding(seq_list):
    '''train noncoding file'''
    
    r_r_counts = [[[0 for i in range(4)] for j in range(4)] for g in range(NUM_STRATIFY)]
        
    for seq in seq_list:
        if STRATIFY:
            gc_content = get_gc_content(seq)
        else:
            gc_content = MIN_GC_CONTENT        
        
        for t in range(len(seq)-1):
            fr = nt_dict.get(seq[t], -1)
            to = nt_dict.get(seq[t+1], -1)
                  
            if fr >= 0 and to >= 0:
                r_r_counts[gc_content - MIN_GC_CONTENT][fr][to] += 1
            
    noncoding_file = open("noncoding", "w")
    
    for gc in range(MIN_GC_CONTENT, MAX_GC_CONTENT+1):
        line = "%s\n" % gc
        noncoding_file.write(line)
        
        if STRATIFY:
            k = gc
        else:
            k = MIN_GC_CONTENT  
        
        for j in range(4):
            total_ct = sum(r_r_counts[k - MIN_GC_CONTENT][j])
            #  print dimer_list[j],
            line = "";
            for i in range(4):
                if total_ct > 0:
                    ct = r_r_counts[k-MIN_GC_CONTENT][j][i]
                    if ct == 0:
                        prob = 0.0001
                    else:
                        prob = round(float(ct) / total_ct, 4)
                        if prob < 0.0001:
                            prob = 0.0001
                else:
                    prob = 0.0001
                line += str(prob)
                line += '\t'
            line = line.strip('\t')
            line += ('\n')
            noncoding_file.write(line)                    
                                  
    noncoding_file.close()
    print "output file produced: noncoding"   
                                 
            
if __name__ == '__main__':
    usage  = "usage: %prog [-i <input sequence file>  | -d <input directory>] [-n <input noncoding file>] [-g]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="input", type = "string", default=None, help="Input gene sequence file.")
    parser.add_option("-d", "--input_dir",  dest="input_dir", type = "string", default=None, help="Directory containing input gene sequence files.")
    parser.add_option("-n", "--noncoding",  dest="noncoding", type = "string", default=None, help="Input noncoding sequence file.")
    parser.add_option("-g", "--gc", dest="gc_content", action="store_true", default=True, help="stratify by gene GC content")
    
    (opts, args) = parser.parse_args()
    
    STRATIFY = opts.gc_content
    msg = "Stratify= %s" % STRATIFY
    
    gene_list = []
    noncoding_list = []
    
    if opts.input:
        if os.path.isfile(opts.input):
            gene_list = parse_input_file(opts.input)    
            if opts.noncoding and os.path.isfile(opts.concoding):
                noncoding_list = parse_input_file(opts.noncoding, noncoding=True)     
        else:
            parser.error("Specified input file not found: " % opts.input)
            
    elif opts.input_dir:
        if os.path.isdir(opts.input_dir):
            os.chdir(opts.input_dir)
            gene_list, noncoding_list = parse_input_dir(opts.input_dir)
        else:
            parser.error("Specified input directory not found: ", opts.input_dir)   
            
    else:
        parser.error("Missing input file or input directory")
        
    if gene_list:
        print "total # of gene sequences=", len(gene_list)
        train_gene_transition_two_way(gene_list)
        train_start_stop_adjacent_prob(gene_list)
    
    if noncoding_list:
        train_non_coding(noncoding_list)
        