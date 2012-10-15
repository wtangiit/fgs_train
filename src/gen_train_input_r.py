#!/usr/bin/env python

'''generate needed input file for FGS training script fgs_train.py. 
  Input: -c <contig> -t <tbl>
  Output:  training data fwd
''' 

import sys, os
from optparse import OptionParser
from Bio import SeqIO

def get_noncoding_seqs(bitmap_dict, direction):
    '''retrieve noncoding area by bitmap'''
    last_bit = 1
    index = 0
    start = 0
    for contig_id in bitmap_dict.keys():
        record_seq = contig_dict[contig_id]
        bitmap = bitmap_dict[contig_id]
        #print len(bitmap), sum(bitmap)
        for bit in bitmap:
            if bit == 1:
                if last_bit == 0:
                    stop = index
                    msg = "%s_%s-%s_%s_noncoding\t%d\t"%(contig_id, start, stop-1, direction, stop-start)
                    
                    if direction == '+':
                        noncoding_seLq = str(record_seq[start:stop])
                    elif direction == '-':
                        noncoding_seq =  str(record_seq[start:stop].reverse_complement())
                    
                    #print "noncoding_seq=", noncoding_seq
                    
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
        
def parse_contigs(contigs_file):
    '''parse contigs file, return a dictionary {congig_id: nt_sequence}'''
    inf = open(contigs_file)
    records=SeqIO.parse(inf, "fasta")
    contig_dict = {}
    
    for seq_record in records:
        #print seq_record
        key = seq_record.id.split("|")[2]
        seq = seq_record.seq
        if not contig_dict.has_key(key):
            contig_dict[key] = seq
        else:
            print "adding a contig with ID %s already existed" % key
    #print "total_contigs:", len(contig_dict.keys())
    inf.close()
    return contig_dict

if __name__ == '__main__':
    usage  = '''usage: %prog -d <input directory> [-c contig_path] [-t tbl_path] [-f] [-r] [-b buffer] \n
    purpose: produces tables of sequence subsets corresponding to genes with upstream and downstream sequences'''
    parser = OptionParser(usage)
    parser.add_option("-c", "--contigs",   dest="contigs", default=None, help="Input contigs file.")
    parser.add_option("-t", "--tbl",     dest="tbl", default=None, help="Input tbl file.")
    parser.add_option("-b", "--buffer",  dest="buf", default=60, help="Forward / reverse buffer region")
    parser.add_option("-d", "--dir",   dest="dir", default=".", help="Input contigs file.")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose [default off]")
    parser.add_option("-f", "--fasta", dest="fasta", action="store_true", default=False, help="Fasta output (default csv)")
    
    (opts, args) = parser.parse_args()
    buf = int(opts.buf)
    
    if not (opts.dir and os.path.isdir(opts.dir)):
        parser.error(("Missing input directory name %s"%(opts.dir, )))
    
    dirname = opts.dir
    
    if opts.contigs==None:
        contig_file = "%s/contigs" % dirname
    else:
        contig_file = contig_file
    if not os.path.isfile(contig_file):
        parser.error("Cannot find input contigs file %s"%(contig_file, ))
        
    if opts.tbl==None:
        tbl_file = "%s/Features/peg/tbl" % dirname
    else:
        tbl_file = tbl_file
    if not os.path.isfile(tbl_file):
        parser.error("Cannot find input tbl file %s"%(tbl_file, ))
        
    print "input config=%s, tbl=%s, buf=%s" % (contig_file, tbl_file, buf)
    upstream = buf
    downstream = buf
    
    if opts.verbose:
        sys.stderr.write("Processing %s and %s... \n" % (contig_file, tbl_file))
    
    contig_dict = parse_contigs(contig_file)
    
    #creating bitmaps marking start stop
    bitmap_fwd_dict = {}
    bitmap_rev_dict = {}
    for key in contig_dict.keys():
        rec_seq = contig_dict[key]
        bitmap_fwd_dict[key] = [0 for i in range(len(rec_seq))]
        bitmap_rev_dict[key] = [0 for i in range(len(rec_seq))]
        
    outfilename_coding = "train_fwd_%s" % dirname
    outfilename_noncoding = "train_noncoding_%s.csv" %  dirname
    
    if opts.fasta:
        outfilename_coding += ".fasta"
    else:
        outfilename_coding += ".csv"
        
    outfile_coding = open(outfilename_coding, "w")
    outfile_noncoding = open("train_noncoding_%s.csv" %  dirname, "w")
    
    tbl_handle = open(tbl_file)
                    
    for line in tbl_handle:
        line = line.strip()
        useful_seg = line.split("|")[-1]
        splits = useful_seg.split("_")
        contig_id = "%s_%s" % (splits[0], splits[1])
        left = int(splits[2])
        right = int(splits[3])
        if left <= right:
            direction = "+"
            start = left
            stop = right
        else:
            direction = '-'
            start = right
            stop = left
            
        record_seq = contig_dict[contig_id]
            
        #print contig_id, direction, start, stop, record_seq[0:10], record_seq[0:10].reverse_complement(), len(record_seq)
        
        if direction == "+":
            if start-1-upstream < 0 or stop+downstream > len(record_seq):
                #sys.stderr.write("problem with %s out of range. \n"%(label))
                continue
            
            msg = ""
            if opts.fasta:
                #print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
                msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(contig_id, start, stop, direction, buf, (stop-start))
            else:
                #print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
                msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(contig_id, start, stop, direction, buf, (stop-start))
            #seq = str(record_seq[(int(start)-1):(int(stop))])
            seq = str(record_seq[(int(start)-1-upstream):(int(stop)+downstream)])
  
            #print seq
            msg += seq
            msg += '\n'
            outfile_coding.write(msg)
            
            #marking coding area
            for i in range(start-1, stop):
                bitmap_fwd_dict[contig_id][i] = 1
                 
        elif direction == "-":
            if start-1-downstream  < 0 or stop + upstream > len(record_seq):
                #sys.stderr.write("problem with %s out of range. \n"%(label))
                pass
            else:
                msg = ""
                if opts.fasta:
                    msg += ">%s_%s-%s_%s_plusminus%d length=%d\n"%(contig_id, start, stop, direction, buf, (stop-start))
                else:
                    msg += "%s_%s-%s_%s_plusminus%d\t%d\t"%(contig_id, start, stop, direction, buf, (stop-start))
                    
                #seq = str(record_seq[int(start)-1:int(stop)].reverse_complement())
                seq = str(record_seq[int(start)-1-downstream:int(stop)+upstream].reverse_complement() )
                
                #print seq
                msg+= seq
                msg += '\n'
                outfile_coding.write(msg)
                
                #marking coding area
                for i in range(start-1, stop):
                    bitmap_rev_dict[contig_id][i] = 1
                    
        else:
            print "invalid direction ", direction
            pass
                
#retrieve noncoding area
    get_noncoding_seqs(bitmap_fwd_dict, '+')
    get_noncoding_seqs(bitmap_rev_dict, '-')
                    
    tbl_handle.close()
    outfile_coding.close()
    outfile_noncoding.close()
    
    if opts.verbose: sys.stderr.write("Done. \n")
