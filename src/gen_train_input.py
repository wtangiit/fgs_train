#!/usr/bin/env python

'''generate needed input file for FGS training script fgs_train.py. 
  Input: <input sequence file>  <input ptt table>
  Output:  train data for training gene and rgene,  training data for training noncoding
''' 

import sys, os
from optparse import OptionParser
from Bio import SeqIO

def get_noncoding_seqs(bitmap, direction, label, record, outfile_noncoding):
    '''retrieve noncoding area by bitmap'''
    last_bit = 1
    index = 0
    start = 0
    for bit in bitmap:
        if bit == 1:
            if last_bit == 0:
                stop = index
                msg = "%s_%s-%s_%s_noncoding\t%d\t"%(label, start, stop-1, direction, stop-start)
                
                if direction == '+':
                    noncoding_seq = str(record.seq[start:stop])
                elif direction == '-':
                    noncoding_seq =  str(record.seq[start:stop].reverse_complement())
                
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
        
def gen_train_seqs(input_fna, tempdata_dir, input_ptt=None, fasta=False, buf=60):
    '''generate csv from raw input fna'''

    
    segs = input_fna.split('/')
    last = segs[-1]
    segs2 = last.split('.')
    seqid = segs2[0]
    
    print "generating .csv for:", seqid
    
    if not input_ptt:
        input_ptt = os.path.splitext(input_fna)[0] + ".ptt"
        print "input_ptt=", input_ptt
        
    if not os.path.isfile(input_ptt):
        print "Missing input ptt file %s"% (input_ptt)
        sys.exit()
        
    print "inputfile=%s, ptt file=%s, buf=%s" % (input_fna, input_ptt, buf)
    upstream = buf
    downstream = buf
    
    in_handle  = open(input_fna)
    ptt_handle = open(input_ptt)
    record=SeqIO.parse(in_handle, "fasta").next()
    
    bitmap_fwd = [0 for i in range(len(record))]
    bitmap_rev = [0 for i in range(len(record))]
    
    outfilename_coding = "%s/train_fwd_%s" % (tempdata_dir, seqid)
    outfilename_noncoding = "%s/train_noncoding_%s.csv" % (tempdata_dir, seqid)
    
    if fasta:
        outfilename_coding += ".fasta"
    else:
        outfilename_coding += ".csv"
        
    outfile_coding = open(outfilename_coding, "w")
    outfile_noncoding = open(outfilename_noncoding, "w")
                    
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
                    #sys.stderr.write("problem with %s out of range. \n"%(label))
                    pass
                else:
                    msg = ""
                    if fasta:
                        #print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
                        msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(label, start, stop, direction, buf, (stop-start))
                    else:
                        #print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
                        msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start))
                    #seq = str(record.seq[(int(start)-1):(int(stop))])
                    seq = str(record.seq[(int(start)-1-upstream):(int(stop)+downstream)])
          
                    #print seq
                    msg += seq
                    msg += '\n'
                    outfile_coding.write(msg)
                    
                    #marking coding area
                    for i in range(start-1, stop):
                        bitmap_fwd[i] = 1
                     
            if direction == "-":
                if start-1-downstream  < 0 or stop + upstream > len(record.seq):
                    #sys.stderr.write("problem with %s out of range. \n"%(label))
                    pass
                else:
                    msg = ""
                    if fasta:
                    #   print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
                        msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(label, start, stop, direction, buf, (stop-start))
                    else:
                    #   print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
                        msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start))
                        
                    #seq = str(record.seq[int(start)-1:int(stop)].reverse_complement())
                    seq = str(record.seq[int(start)-1-downstream:int(stop)+upstream].reverse_complement() )
                    
                    #print seq
                    msg+= seq
                    msg += '\n'
                    outfile_coding.write(msg)
                    
                    #marking coding area
                    for i in range(start-1, stop):
                        bitmap_rev[i] = 1    
                    
    #retrieve noncoding area   ???needcheck
    get_noncoding_seqs(bitmap_fwd, '+', label, record, outfile_noncoding)
    get_noncoding_seqs(bitmap_rev, '-', label, record, outfile_noncoding)
                    
    in_handle.close()
    outfile_coding.close()
    outfile_noncoding.close()
    

if __name__ == '__main__':
    usage  = '''usage: %prog -i <input sequence file> -p <input ptt> [-f] [-r] [-b buffer] \n
    purpose: produces tables of sequence subsets corresponding to genes with upstream and downstream sequences'''
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",   dest="input", default=None, help="Input sequence file.")
    parser.add_option("-p", "--ptt",     dest="ptt", default=None, help="Input ptt table.")
    parser.add_option("-b", "--buffer",  dest="buf", default=60, help="Forward / reverse buffer region")
    parser.add_option("-f", "--fasta", dest="fasta", action="store_true", default=False, help="Fasta output (default csv)")
    parser.add_option("-t", "--tab", dest="tabsep", action="store_true", default=False,  help="Tab separated")
    
    (opts, args) = parser.parse_args()
    buf = int(opts.buf)
    if not (opts.input and os.path.isfile(opts.input) ):
        parser.error("Missing input file %s"%(opts.input, ))
        sys.exit(1)
            
    gen_train_seqs(opts.input)
    
         
  
    