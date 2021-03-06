#!/usr/bin/env python

'''generate needed input file for FGS training script fgs_train.py. 
  Input: -c <contig> -t <tbl>
  Output:  training data fwd
''' 

import sys, os
from optparse import OptionParser
from Bio import SeqIO
import fgs_train

def get_noncoding_seqs(bitmap_dict, direction, contig_dict, outfile_noncoding):
    '''retrieve noncoding area by bitmap'''

    noncoding_seqs = []
    
    for contig_id in bitmap_dict.keys():
        last_bit = 1
        index = 0
        start = 0
        record_seq = contig_dict[contig_id]
        bitmap = bitmap_dict[contig_id]
        #print len(bitmap), sum(bitmap)
        for bit in bitmap:
            if bit == 1:
                if last_bit == 0:
                    stop = index
                    msg = "%s_%s-%s_%s_noncoding\t%d\t"%(contig_id, start, stop-1, direction, stop-start)
                    
                    if direction == '+':
                        noncoding_seq = str(record_seq[start:stop])
                    elif direction == '-':
                        noncoding_seq =  str(record_seq[start:stop].reverse_complement())
                    msg += noncoding_seq
                    msg += '\n'
                    if opts.writecsv:
                        outfile_noncoding.write(msg)

                    noncoding_seqs.append(noncoding_seq)
            elif bit == 0:
                if last_bit == 1:
                    start = index
                pass
            index += 1
            last_bit = bit
    return noncoding_seqs
        
def parse_contigs(contigs_file):
    '''parse contigs file, return a dictionary {congig_id: nt_sequence}'''
    inf = open(contigs_file)
    records=SeqIO.parse(inf, "fasta")
    contig_dict = {}
    
    for seq_record in records:
        key = seq_record.id
        seq = seq_record.seq
        if not contig_dict.has_key(key):
            contig_dict[key] = seq
        else:
            print "adding a contig with ID %s already existed" % key
    print "total_contigs:", len(contig_dict.keys())
    inf.close()
    return contig_dict

def gen_train_seqs(dirname, contig_file=None, tbl_file=None):
    coding_seqs = []
    noncoding_seqs = []
    
    if contig_file==None:
        contig_file = "%s/contigs" % dirname
    
    if not os.path.isfile(contig_file):
        print "Cannot find input contigs file %s"%(contig_file, )
        sys.exit(1)
        
    if tbl_file==None:
        tbl_file = "%s/Features/peg/tbl" % dirname
    if not os.path.isfile(tbl_file):
        print "Cannot find input tbl file %s"%(tbl_file,)
        sys.exit(1)
        
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
    
    if opts.writecsv:
        outfile_coding = open(outfilename_coding, "w")
        outfile_noncoding = open("train_noncoding_%s.csv" %  dirname, "w")
    else:
        outfile_noncoding = None
    
    tbl_handle = open(tbl_file)
                    
    for line in tbl_handle:
        line = line.strip()
        useful_seg = line.split()[-1]
        splits = useful_seg.split("_")
        contig_id = "_".join(splits[0:-2])
        #print contig_id
        left = int(splits[-2])
        right = int(splits[-1])
        if left <= right:
            direction = "+"
            start = left
            stop = right
        else:
            direction = '-'
            start = right
            stop = left
        
        if contig_dict.has_key(contig_id):
            record_seq = contig_dict[contig_id]
        else:
            #print "warning: table entry not found in contig file:", line
            continue
            
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
            if opts.writecsv:
                outfile_coding.write(msg)
            coding_seqs.append(seq)
            
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
                if opts.writecsv:
                    outfile_coding.write(msg)
                coding_seqs.append(seq)
                
                #marking coding area
                for i in range(start-1, stop):
                    bitmap_rev_dict[contig_id][i] = 1
                    
        else:
            print "invalid direction ", direction
            pass
                
#retrieve noncoding area
    noncoding_seqs.extend(get_noncoding_seqs(bitmap_fwd_dict, '+', contig_dict, outfile_noncoding))
    noncoding_seqs.extend(get_noncoding_seqs(bitmap_rev_dict, '-', contig_dict, outfile_noncoding))
                    
    tbl_handle.close()
    
    if opts.writecsv:
        outfile_coding.close()
        outfile_noncoding.close()
    
    return coding_seqs, noncoding_seqs

if __name__ == '__main__':
    usage  = '''usage: %prog -d <input directory> [-c contig_path] [-t tbl_path] [-f] [-r] [-b buffer] \n
    purpose: produces tables of sequence subsets corresponding to genes with upstream and downstream sequences'''
    parser = OptionParser(usage)
    parser.add_option("-c", "--contigs",   dest="contigs", default=None, help="Input contigs file.")
    parser.add_option("-t", "--tbl",     dest="tbl", default=None, help="Input tbl file.")
    parser.add_option("-b", "--buffer",  dest="buf", default=60, help="Forward / reverse buffer region")
    parser.add_option("-d", "--dir",   dest="dir", default=".", help="Input contigs file.")
    parser.add_option("-w", "--writecsv", dest="writecsv", action="store_true", default=False, help="write csv file [default off]")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose [default off]")
    parser.add_option("-f", "--fasta", dest="fasta", action="store_true", default=False, help="Fasta output (default csv)")
    parser.add_option("-p", "--prefix", dest="prefix", default="", help="prefix string for output")
    
    (opts, args) = parser.parse_args()
    buf = int(opts.buf)
    
    if not (opts.dir and os.path.isdir(opts.dir)):
        parser.error(("Missing input directory name %s"%(opts.dir, )))
    
    dirname = opts.dir
    
    if opts.contigs==None:
        contig_file = "%s/contigs" % dirname
    else:
        contig_file = opts.contigs
    if not os.path.isfile(contig_file):
        parser.error("Cannot find input contigs file %s"%(contig_file, ))

    if opts.tbl==None:
        tbl_file = "%s/Features/peg/tbl" % dirname
    else:
        tbl_file = opts.tbl
    if not os.path.isfile(tbl_file):
        parser.error("Cannot find input tbl file %s"%(tbl_file, ))
        
    print "input config=%s, tbl=%s, buf=%s" % (contig_file, tbl_file, buf)
    upstream = buf
    downstream = buf
    
    if opts.verbose:
        sys.stderr.write("Processing %s and %s... \n" % (contig_file, tbl_file))
        
    coding_seqs, noncoding_seqs = gen_train_seqs(dirname, contig_file, tbl_file)
    
    for seq in noncoding_seqs:
#        print len(noncoding_seqs)
        pass
            
    fgs_train.train_gene_transition_two_way(coding_seqs, opts.prefix)
    fgs_train.train_start_stop_adjacent_prob(coding_seqs, opts.prefix)
    fgs_train.train_non_coding(noncoding_seqs, opts.prefix)
    
    if opts.verbose: sys.stderr.write("Done. \n")
