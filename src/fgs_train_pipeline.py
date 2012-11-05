#!/usr/bin/env python
'''FGS training pipeline
This is a wrapper script that invokes other fgs training scripts as a pipeline. 
by invoking this single script only we can get a training output data from the 
proper prepared input data. 
''' 
import os
import sys
from optparse import OptionParser

import ConfigParser
 
import subprocess
 
import gen_train_input
import fgs_train 
   
def parse_genome_list(list_file):
  
    '''parse the genome list file, return a python list of genome id'''
 
    infile = open(list_file, "r")
    id_list = []
    for line in infile:
        segs = line.split()
        id_list.append(segs[0].split('.')[0])
  
    return id_list
   
def iterate_director(dir_path, key_list=None):

    '''iterate all files in the specified directory or files specified in key_list'''
 
    fna_file_list = []
    if key_list:
        for key in key_list:
            fna_file = key + ".fna"
            fna_file_path = dir_path + '/'  + fna_file
            if os.path.isfile(fna_file_path):
                fna_file_list.append(fna_file_path)
                print fna_file_path
    else:
        print "iterate all: ", dir_path
        listing = os.listdir(dir_path)
        for infile in listing:
            if infile.split('.')[-1] == "fna":
                fna_file_path = dir_path + '/' + infile
                fna_file_list.append(fna_file_path)
    return fna_file_list

if __name__ == '__main__':
    usage  = '''usage: %prog -d <dir> -l <list_file>'''
    parser = OptionParser(usage)
    parser.add_option("-d", "--dir", dest="dir", default=None, help="recompile the specified tag")
    parser.add_option("-l", "--list",   dest="list", default=None, help="config file")
    
    (opts, args) = parser.parse_args()
    
    if not (opts.dir and os.path.isdir(opts.dir)):
        parser.error("Missing directory path: %s" % (opts.dir, ))
        sys.exit(1)
        
    data_dir = opts.dir.rstrip('/')
    print "data_dir=", data_dir
        
    train_name = "ALL_" + opts.dir.split('/')[-1]
    
    genome_list = []
    
    if opts.list:
        if not os.path.isfile(opts.list):
            parser.error("Missing list file: %s" % (opts.list, ))
            sys.exit(1)
        list_file = opts.list
        print "list_file=", list_file
        genome_list = parse_genome_list(list_file)
        train_name = list_file.split('/')[-1].split('.')[0]
    else:
        print "list file not specified, will iterate all sequence files within input directory"
    
    #parent_dir = os.path.dirname(data_dir)
    tempdata_dir = "fgs_train_output_" + train_name
    print "fgs_train_output dir: ", tempdata_dir
    
    if not os.path.exists(tempdata_dir):
        os.makedirs(tempdata_dir)
        
    if len(genome_list) > 0:
        fna_file_list = iterate_director(data_dir, genome_list)
    else:
        fna_file_list = iterate_director(data_dir)
    
    for fna_file in fna_file_list:
        gen_train_input.gen_train_seqs(fna_file, tempdata_dir)
        
    fgs_train.fgs_train_dir(tempdata_dir)
