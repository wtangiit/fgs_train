#!/usr/bin/env python
import os
from optparse import OptionParser

def parse_input_file(filename):
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
                if seq_len < 123:
                    print "**********************************************"
                    print "Warning!!:the input data contains invalid data, line %d, sequence length=%s" % (linenumber, seq_len)
                    print "The invalid sequence is thrown out to continue, but replacing input data and re-training is suggested." 
                    print "**********************************************"
                    continue
                seq_lists.append(splits[2])
    return seq_lists


if __name__ == '__main__':
    usage  = '''usage: %prog -i <input sequence file>'''
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",   dest="input", default=None, help="Input sequence file.")
    
    (opts, args) = parser.parse_args()
    
    if not (opts.input and os.path.isfile(opts.input) ):
        parser.error("Missing input file %s"%(opts.input, ))
        
    seq_lists = parse_input_file(opts.input)
    
    print "total number of sequences:", len(seq_lists)
    
    start_dict = {}
    stop_dict = {}
    
    for seq in seq_lists:
        start_codon = seq[60:63]
        stop_codon = seq[-63:-60]
        
        if start_codon in start_dict.keys():
            start_dict[start_codon] += 1
        else:
            start_dict[start_codon] = 1
            
        if stop_codon in stop_dict.keys():
            stop_dict[stop_codon] += 1
        else:
            stop_dict[stop_codon] = 1
            
    print "start codon statistics:"
    for codon in start_dict.keys():
        print codon, ":", start_dict[codon] 
        
    print "stop codon statistics:"
    for codon in stop_dict.keys():
        print codon, ":", stop_dict[codon]