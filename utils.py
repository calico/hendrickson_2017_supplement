#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:43:43 2017

@author: ilyasoifer
"""

import numpy as np
import scipy
import types
import pickle, itertools, sys
from pysam import AlignmentFile
#import pysam
import varutils, genomic_utils

def check_good_read(res):
    '''
    Checks if the read is 
    a. Have a high mapping quality
    b. Mate pair is mapped
    c. Forward-Reverse orientation
    
    Input: 
        res - pysam Alignment
    Output: 
        boolean True/False
    '''
    
    if res.is_secondary:
        return False
    if res.is_unmapped or res.mapping_quality < 30 : 
        return False
    if res.is_paired : 
            if res.is_read2:  # working only on read1 
                return False
            if res.is_unmapped or res.mate_is_unmapped:
                return False
            if not res.is_reverse ^ res.mate_is_reverse : 
                return False
            if not check_fr(res):
                return False
    return True
    
def check_fr( res ):
    '''
    Checks if the alignment is Forward-Reverse
    '''
    
    if not res.is_reverse and res.reference_start <= res.next_reference_start:
        return True
    elif res.is_reverse and res.reference_start >= res.next_reference_start : 
        return True
    return False


    
def pesam2bed( input_file, chrlen_file, output_file, logfile, errfile,
              filter_function = check_good_read):
    '''
    This 
    '''
    
    samfile = AlignmentFile(input_file)

    # check good_read also checks if this is R1... 
    gsamfile = itertools.groupby( samfile, lambda x : x.query_name )
    take_read = [ y for x in itertools.imap(lambda x : filter_function(x[1]),
                                                gsamfile ) for y in x ]
    samfile.reset()
   
    logfile.write("%d out of %d reads passed filtering\n"%(sum(take_read), len(take_read)/2))
    logfile.flush()
    chrlens = genomic_utils.get_chr_lens(chrlen_file)
    
    outfile = open(output_file,'w')
    
    count = 0
    for i,r in enumerate(samfile):
        count += 1
        if count % 1000000 == 0:
            logfile.write("pesam2bed: %s, Read # %d\n"%(input_file,count))
            logfile.flush()

        read_1 = r
        if take_read[i] == 0 :
            continue
        

        if not read_1.is_reverse :
            rstart = read_1.reference_start

            start = max( 0, rstart )        
            rend = rstart + read_1.template_length

            end = min( rend, chrlens[ samfile.getrname(read_1.reference_id) ] )
            if start > end :
                continue
            try:
                outfile.write(samfile.getrname(read_1.reference_id))
                outfile.write('\t')
                outfile.write('%d\t%d\t%s\t%d\t+\n'%(start,end, read_1.qname,
                                                     take_read[i]*100) )
            except ValueError :
                errfile.write("***** pesam2bed failed with the following value\n")
                errfile.write(samfile.getrname(read_1.reference_id))
                errfile.write(' %d\t%d\t%s\t%d\t-\n'%(start,end, read_1.qname,
                                                      take_read[i]*100))
                errfile.flush()
                sys.exit(-1)
            
            
        else:
            rstart = read_1.reference_end - abs(read_1.template_length)

            start = max(0, rstart)
            rend = read_1.reference_end
            end = min( rend, chrlens[ samfile.getrname(read_1.reference_id) ])
            if start > end:
                continue
            try:
                outfile.write(samfile.getrname(read_1.reference_id))
                outfile.write('\t')
                outfile.write('%d\t%d\t%s\t%d\t-\n'% (start,end, read_1.qname,
                                                      take_read[i]*100) )
            except ValueError :
                errfile.write("***** pesam2bed failed with the following value\n")
                errfile.write(samfile.getrname(read_1.reference_id))
                errfile.write(' %d\t%d\t%s\t%d\t-\n'%(start,end, read_1.qname,
                                                      take_read[i]*100))
                errfile.flush()
                sys.exit(-1)

    outfile.close()




def get_chr_coverage( fragments, size_range, chrlen_file,
                     out=sys.stdout ): 
    '''
    
    '''
    chrlens = genomic_utils.get_chr_lens(chrlen_file)
    chr_coverage = {}
    for chrom in chrlens:
        chr_coverage[chrom] = scipy.sparse.dok_matrix((1,chrlens[chrom]))

    with open(fragments) as fragments_file :
        it1,it2 = itertools.tee(genomic_utils.parse_bed_line_gen(fragments_file))
        # 4 is the shift between the insert point to the actual site
        itstart = itertools.ifilter(lambda x:x.length >= size_range[0] \
                                and x.length <= size_range[1], \
                                itertools.imap(lambda x: genomic_utils.Position((x[0],x[1] \
                                + 4,x[2]-x[1]-8,x[3]), weight=x[5]),it1 ))
        # BED is half-open, so the end does not belong to the fragment
        itend    = itertools.ifilter(lambda x:x.length >= size_range[0] \
                                and x.length <= size_range[1], \
                                itertools.imap(lambda x: genomic_utils.Position((x[0],x[2]-5,\
                                x[2]-x[1]-8,x[3]), weight=x[5]),it2 )) 
        
        for cc,pos in enumerate( varutils.roundrobin( itstart, itend ) ):

            if cc%500000==0 :
                print >>out, "get_chr_coverage:", fragments, "position", pos,cc
            chr_coverage[ pos.chromosome ][ 0,pos.pos ]+=pos.weight
            
    for chrom in chrlens:
        chr_coverage[ chrom ] = chr_coverage[chrom].tocsr()
            
    return chr_coverage


def get_ins_coverage( fragments_bed, intervals, size_range = None, \
                     chrlen_file = None, norm_factor = None, \
                     blacklist = '', normlist = '', pseudocount = 0): 
    '''
    Input: 
        fragments  
        intervals 
        size_range
        chrlen_file
        norm_factor
        blacklist
        normlist
        pseudocount
    if  size_range a string ('low','high' or 'all'), 
    fragments_bed should be .pickle file that the pipeline generates
    Otherwise it is assumed to be a bed file of 
    aligned reads and chr coverage is calculated anew 
    
    '''
    pseudocount = float(pseudocount)
    if type(fragments_bed)!=types.ListType:
        fragments_bed = [ fragments_bed ]
    if type(size_range) == type('low') and size_range=='low':
        chr_coverages = [ pickle.load(open(x))['low'] for x in fragments_bed ]
    elif type(size_range) == type('high') and size_range == 'high':   
        chr_coverages =  [ pickle.load(open(x))['high'] for x in fragments_bed ]
    elif type(size_range) == type('all') and size_range == 'all':   
        chr_coverages =  [ pickle.load(open(x))['all'] for x in fragments_bed ]

    elif type(fragments_bed[0]) == type({}):
        chr_coverages = fragments_bed
    else:
        chr_coverages = [ get_chr_coverage( x, size_range, chrlen_file ) 
                                for x in fragments_bed ]
    
    chr_coverage = unite_chr_coverage( chr_coverages )
    if not norm_factor:
        if blacklist: 
            norm_factor = get_norm_factor( chr_coverage, blacklist=open(blacklist ))
        elif normlist:
            norm_factor = get_norm_factor( chr_coverage, whitelist=open(normlist ))
        else:
            norm_factor = get_norm_factor( chr_coverage,  '', '')
    
    intervals = map(lambda x: genomic_utils.Interval(x,True), intervals)
    interval_dct = {}
    count = 0
    intervals  = sorted(intervals, key = lambda x: x.chromosome )
    gintervals = itertools.groupby(intervals, lambda x:x.chromosome)
    for grp in gintervals :
        if grp[0] not in chr_coverage :
                continue
        if type( chr_coverage[grp[0]] ) == type( np.ndarray((0,0))):
            chromosome = chr_coverage[grp[0]]
        else: 
            chromosome = chr_coverage[grp[0]].toarray()
        if chromosome.ndim == 1 :
            chromosome = chromosome.reshape(1,-1)

        for interval in grp[1]:
            count += 1
            if interval.start < 0 or interval.end > chromosome.shape[1] :
                continue
            
            new_data = chromosome[ 0,interval.start:interval.end ].squeeze().reshape(-1)
            if interval.strand == '+':
                if hasattr(interval,'gene'):
                    interval_dct[interval.gene] = \
                                (new_data + pseudocount/len(interval)) * norm_factor
                else:
                    interval_dct[hash(interval)] = \
                                 ( new_data + pseudocount/len( interval ) )  * norm_factor
                    
            else:
                if hasattr(interval, 'gene'):
                    interval_dct[interval.gene] = \
                                ( new_data[::-1] + pseudocount/len(interval) ) \
                                * norm_factor
                else: 
                    interval_dct[hash(interval)] = \
                                 ( new_data[::-1] + pseudocount/len(interval) )\
                                 * norm_factor
    return interval_dct
    

def check_good_read(res):
    if res.is_secondary:
        return False
    if res.is_unmapped or res.mapping_quality < 30 : 
        return False
    if res.is_paired : 
            if res.is_read2:  # working only on read1 
                return False
            if res.is_unmapped or res.mate_is_unmapped:
                return False
            if not res.is_reverse ^ res.mate_is_reverse : 
                return False
            if not check_fr(res):
                return False
    return True
    





def unite_chr_coverage(coverages):
    if not coverages:
        return None
    result = coverages[0]
    for c in coverages[1:]:
        for key in result.keys():
            result[key]+=c[key]
    return result

    
def get_norm_factor( chrcov, blacklist = [], whitelist = [] ):
    '''
    Get the normalization factor 
    
    Input: 
        chrcov - dictionary of nucleotide coverages for each chromosome
        
        blacklist - optional list of lines from blacklist file, or handle to blacklist file 
        lines are tab-separated: chromosome start end (optional fields after are ignored ),
        start, end are 1-based, interval is closed
        optionally line could be chromosome (without start and end)
        
        whitelist - similar to blacklist but the normalization coverage will be calculated 
        only based on intervals in the whitelist
    Output: 
        
        normalization factor, that the coverages should be multiplied by so that the 
        sum of all coverages outside blacklist is 100 ( or sum of all coverages in whitelist are 100)
            
    '''
    vals = [ chrcov[x].sum() for x in chrcov ]
    if whitelist:
        whitelist_flag = True
    else: 
        whitelist_flag = False
        
    black_lsp = [ x.strip().split("\t") for x in blacklist ]
    white_lsp = [ x.strip().split("\t") for x in whitelist]
    if not black_lsp:
        black_lsp = white_lsp
        
    black_cvg = 0 
    for line in black_lsp :
        chromosome = line[0]
        if len( line ) > 1 :
            start = int(line[1])-1
            end = int(line[2])+1
        else:
            start = 0 
            end = chrcov[chromosome].shape[1]
        black_cvg += (chrcov[chromosome][0,start:end]).sum()
    if whitelist_flag : 
        return 100./(black_cvg)
    return 100./( sum(vals) - black_cvg )

