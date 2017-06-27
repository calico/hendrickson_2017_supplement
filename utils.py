#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import scipy.sparse
import pickle, itertools, sys, os
from pysam import AlignmentFile
import varutils, genomic_utils
def check_good_read(res):
    '''
    Checks if the read is 
    a. Have a high mapping quality
    b. Mate pair is mapped
    c. Forward-Reverse orientation
    
    Input: 
        res - pysam.AlignedSegment
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
    Input: 
        res - pysam.AlignedSegment
    '''
    
    if not res.is_reverse and res.reference_start <= res.next_reference_start:
        return True
    elif res.is_reverse and res.reference_start >= res.next_reference_start : 
        return True
    return False


    
def pesam2bed( input_file, chrlen_file, output_file, logfile, errfile,
              filter_function = check_good_read):
    '''
    Converts aligned BAM file into BED file of fragments filtering reads that do not pass the filtering criteria
    
    Input: 
        input_file - aligned (unsorted) paired end BAM
        chrlen_file - TSV file of format chromosome<TAB>length
        output_file - output file (BED file with start and end corresponding to both ends of the PE fragment)
        logfile, errfile - logging files
        filter_function - function that accepts an aligned read (pysam.AlignedSegment) and returns true, false or weight of the fragment
                            Default: utils.check_good_reads
                            Note: the function can return weight, if trying to accommodate reads from repetitive sequences
    
    Output: 
        None (written into output file)
    
    See also: 
        check_good_read
        
    '''
    
    samfile = AlignmentFile(input_file)

    take_read = [ x for x in itertools.imap(filter_function, samfile)]
    # check good_read also checks if this is R1... 
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
    Reads PEBED with ATACseq insertions and generates a dictionary of per-nucleotide ATAC insertion density. 
    Ends of fragments are shifted by 4 and 5 nucleotides to accommodate for the Tn5 binding site location. 
    
    Input: 
        fragments - PEBED file
        size_range - limits for fragment sizes to count as valid insertions [min_length, max_length]
        chrlen_file - TSV file of format chromosome<TAB>length
        out - handle to write the log\err to
    
    Output: 
        Dictionary (with keys=chromosomes) of vectors with number of insertions at every nucleotide
    '''
    
    chrlens = genomic_utils.get_chr_lens(chrlen_file)
    chr_coverage = {}
    for chrom in chrlens:
        chr_coverage[chrom] = scipy.sparse.dok_matrix((1,chrlens[chrom]))

    with open(fragments) as fragments_file :
        it1,it2 = itertools.tee(genomic_utils.parse_bed_line_gen(fragments_file))
        itstart = itertools.ifilter(lambda x:x.length >= size_range[0] \
                                and x.length <= size_range[1], \
                                itertools.imap(lambda x: genomic_utils.Position((x[0],x[1] \
                                + 4,x[2]-x[1]-8,x[3]), weight=x[5]),it1 ))
        itend    = itertools.ifilter(lambda x:x.length >= size_range[0] \
                                and x.length <= size_range[1], \
                                itertools.imap(lambda x: genomic_utils.Position((x[0],x[2]-5,\
                                x[2]-x[1]-10,x[3]), weight=x[5]),it2 )) 
        
        for cc,pos in enumerate( varutils.roundrobin( itstart, itend ) ):
    
            if cc%500000==0 :
                print >>out, "get_chr_coverage:", fragments, "position", pos,cc
            chr_coverage[ pos.chromosome ][ 0,pos.pos ]+=pos.weight
            
    for chrom in chrlens:
        chr_coverage[ chrom ] = chr_coverage[chrom].tocsr()
            
    return chr_coverage

def get_chr_coverage_vector(chrcov, smoothing_window=1,
                            smoothing_profile='hanning', 
                            logfile=open(os.devnull,'w')):
    '''
    Converts the vector of insert densities to list of positions with non-zero insert density, potentially after smoothing.
    This is useful for writing the data into a file in bedGraph/bigWig formats. 
    
    See also: 
        pipeline.generate_tracks
        
    Input: 
        chrcov - chromosome coverage dictionary (output of get_chr_coverage)
        smoothing_window - size of smoothing window (default: 1, no smoothing)
        smoothing_profile - profile of smoothing window (flat, hanning, hamming, bartlett, blackman)
        logfile - handle of the log file (default - devnull, no logging)
    '''
    
    res = {}
    for x in chrcov:
        if smoothing_window > 1:
            logfile.write("Smoothing %s\n"%x)
            chrc = chrcov[x].toarray()
            chrc = varutils.smooth(chrc[0,:], smoothing_window, window=smoothing_profile)
            nz = chrc.nonzero()
            res[x] = (nz[0], chrc[nz])
        else:
            res[x] = scipy.sparse.find(chrcov[x])[1:]
       
    return res



def get_ins_coverage( fragments_file, intervals, norm_factor = None, \
                     blacklist = '', whitelist = '', pseudocount = 0): 
    '''
    Gets a set of intervals and calculates the depth normalized insert coverage (i.e. insertion density) on each interval
    The depth normalization is done so that either the total number of insertions in the whitelist intervals is constant (if whitelist is given)
    or the total number of insertions everywhere outside the blacklist is constant. 
    
    Input: 
        
        fragments - python pickle file (dictionary with vectors of insert numbers on each chromosome), output of 
                    utils.get_chr_coverage
        intervals - list of genomic_utils.Intervals to calculate insertion density on. 
        norm_factor - normalization factor to convert the coverage to insertion density 
                    (optional, would be calculated internally using blacklist or whitelist if not given)
        blacklist - optional blacklist file: normalization would be calculated on all the genome except for 
                    blacklisted intervals.
        whitelist - similar to blacklist but the normalization coverage will be calculated 
                    only based on intervals in the whitelist

        pseudocount - pseudocount is optioanally added to each 

        Note:  In the whitelist/blacklist lines are tab-separated: chromosome start end (optional fields after are ignored ),
                    start, end are 1-based, interval is closed
                    optionally line could be chromosome (without start and end)

    Output: 

        Dictionary of normalized insertion coverage vectors. Keys are either "gene" property of intervals in the list 
        or hash of the interval if they have no "gene" property. 
    
    See also:     utils.get_norm_factor, pipeline.insert_coverage_per_gene
        
    '''
    
    pseudocount = float( pseudocount )
    
    chr_coverage = pickle.load( open( fragments_file ) )
    
    if not norm_factor:
        if blacklist: 
            norm_factor = get_norm_factor( chr_coverage, blacklist=open(blacklist ))
        elif whitelist:
            norm_factor = get_norm_factor( chr_coverage, normlist=open(whitelist))
        else:
            norm_factor = get_norm_factor( chr_coverage,  '', '')
    
    intervals = map( lambda x: genomic_utils.Interval(x,True), intervals )
    interval_dct = {}
    count = 0
    intervals  = sorted( intervals, key = lambda x: x.chromosome )
    gintervals = itertools.groupby( intervals, lambda x:x.chromosome )

    for grp in gintervals :
        if grp[0] not in chr_coverage :
                continue
            
        chromosome = chr_coverage[grp[0]].toarray()

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
    
    
def write_chr_coverage_vectors( chrcov_vector, output_file,mode='dat' ):
    '''
    Write chromosome coverage either as dat or bedGraph format
    Input: 
        chrcov_vector - output of utils.get_chr_coverage_vector
        output_file - name of output file
        mode - bedGraph or dat
    '''
    out = open(output_file,'w')
    chr_sorted = sorted(chrcov_vector.keys())
    for x in chr_sorted:
        for y in range(len(chrcov_vector[x][0])):
            if mode == 'dat':
                out.write('%s\t%d\t%d\n'%(x, chrcov_vector[x][0][y], chrcov_vector[x][1][y]))
            elif mode == 'bedGraph': 
                out.write('%s\t%d\t%d\t%G\n'%(x, chrcov_vector[x][0][y], chrcov_vector[x][0][y]+1,chrcov_vector[x][1][y]))
                
    out.close()


def get_norm_factor( chrcov, blacklist = [], whitelist = [] ):
    '''
    Get the normalization factor that should multiply the insert density coverages so that 
    the coverages on whitelist or outside of blacklist sum up to 100.
    
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

