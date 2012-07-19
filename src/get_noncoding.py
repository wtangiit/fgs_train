#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from Bio import SeqIO

def get_noncoding_seqs(bitmap, direction):
    '''retrieve noncoding area by bitmap'''
    last_bit = 1
    index = 0
    start = 0
    for bit in bitmap:
        if bit == 1:
            if last_bit == 0:
                stop = index
                msg = "%s_%s-%s_%s_noncoding\t%d\t"%(label, start, stop-1, direction, stop-start)
                noncoding_seq = str(record.seq[start:stop])
                msg += noncoding_seq
                msg += '\n'
                outfile_noncoding.write(msg)
            pass
        
        elif bit == 0:
            if last_bit == 1:
                start = index
            pass
        index += 1
        last_bit = bit       
        

if __name__ == '__main__':
    usage  = '''usage: %prog -i <input sequence file> -p <input ptt> [-f] [-r] [-b buffer] \n
    purpose: produces tables of sequence subsets corresponding to genes with upstream and downstream sequences'''
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",   dest="input", default=None, help="Input sequence file.")
    parser.add_option("-p", "--ptt",     dest="ptt", default=None, help="Input ptt table.")
    parser.add_option("-b", "--buffer",  dest="buf", default=60, help="Forward / reverse buffer region")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
    parser.add_option("-f", "--fasta", dest="fasta", action="store_true", default=False, help="Fasta output (default csv)")
    parser.add_option("-t", "--tab", dest="tabsep", action="store_true", default=False,  help="Tab separated")
    
    (opts, args) = parser.parse_args()
    buf = int(opts.buf)
    if not (opts.input and os.path.isfile(opts.input) ):
        parser.error("Missing input file %s"%(opts.input, ))
    if not (opts.ptt and os.path.isfile(opts.ptt) ):
        parser.error("Missing input file %s"%(opts.ptt, ))
    upstream = buf
    downstream = buf
    if opts.verbose: 
        sys.stderr.write("Processing %s and %s... \n"%(opts.input, opts.ptt))
    in_handle  = open(opts.input)
    ptt_handle = open(opts.ptt)
    record=SeqIO.parse(in_handle, "fasta").next()
    
    bitmap_fwd = [0 for i in range(len(record))]
    bitmap_rev = [0 for i in range(len(record))]
    
    outfilename_coding = "train_fwd_%s" % opts.input
    outfilename_noncoding = "train_noncoding_%s.csv" %  opts.input
    
    if opts.fasta:
        outfilename_coding += ".fasta"
    else:
        outfilename_coding += ".csv"
        
    outfile_coding = open(outfilename_coding, "w")
    outfile_noncoding = open("train_noncoding_%s.csv" %  opts.input, "w")
                    
    for line in ptt_handle:
        fields = line.split()
        field1 = fields[0]
        field2 = field1.split(".")
        if len(field2) == 3 :
            start = int(field2[0])
            stop  = int(field2[2])
            direction = fields[1]
            try:
                label= fields[5]
            except:
                label=""
            if direction == "+":
                if start-1-upstream < 0 or stop+downstream > len(record.seq):
                    sys.stderr.write("problem with %s out of range. \n"%(label))
                else:
                    msg = ""
                    if opts.fasta:
                        #print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
                        msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(label, start, stop, direction, buf, (stop-start))
                    else:
                        #print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
                        msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start))
                    seq = str(record.seq[(int(start)-1):(int(stop))])
          
                    #print seq
                    msg += seq
                    msg += '\n'
                    outfile_coding.write(msg)
                    
                    #marking coding area
                    for i in range(start-1, stop):
                        bitmap_fwd[i] = 1
                     
            if direction == "-":
                if start-1-downstream  < 0 or stop + upstream > len(record.seq):
                    sys.stderr.write("problem with %s out of range. \n"%(label))
                else:
                    msg = ""
                    if opts.fasta:
                    #   print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
                        msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(label, start, stop, direction, buf, (stop-start))
                    else:
                    #   print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
                        msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start))
                        
                    seq = str(record.seq[int(start)-1:int(stop)].reverse_complement())
                    
                    #print seq
                    msg+= seq
                    msg += '\n'
                    outfile_coding.write(msg)
                    
                    #marking coding area
                    for i in range(start-1, stop):
                        bitmap_rev[i] = 1    
                    
    #retrieve noncoding area
    get_noncoding_seqs(bitmap_fwd, '+')
    get_noncoding_seqs(bitmap_rev, '-')
                    
    in_handle.close()
    outfile_coding.close()
    outfile_noncoding.close()
    
    if opts.verbose: sys.stderr.write("Done. \n")
