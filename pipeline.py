#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import utils, genomic_utils
import pickle, tempfile, os, subprocess, sys
import pandas as pd


def bam2intervals( input_file, output_file, chrlen_file ):
    '''
    Converts aligned reads (BAM file) into PEBED file after filtering reads with 
    bad alignments. 
    
    See also: 
        utils.pesam2bed
    
    Input: 
        input_file: Unsorted BAM file with alignments of paired end reads
        output_file - three file names with output BED file, logfile and errfile names
        chrlen_file - TSV file with columns: chromosome, chromosome length
    '''
    
    input_file = input_file
    outfile, logfile, errfile = output_file
    logfile = open(logfile,'w')
    errfile = open(errfile,'w')
    
    utils.pesam2bed( input_file, chrlen_file,
                                outfile, logfile, errfile)
    logfile.close()
    errfile.close()
    

def get_coverage_vectors_atac_seq( input_file, output_file, chrlens ):
    
    '''
    get_coverage_vectors_atac_seq - gets the coverage of ATACseq inserts
    
    Input:
        input_file - bed file of inserts (output of bam2intervals).
        output_file[0] - pickle file with dictionary (per-chromosome) of insert coverages.
        output_file[1] - log file.
        chrlen_file - TSV file with columns: chromosome, chromosome length.
    
    Output: 
        Returns insertion coverage per nucleotide (in a dictionary with a key for each chromosome).
    See also: 
        utils.get_chr_coverage
        
    '''
    
    with open( output_file[1], 'w' ) as log :
        fragments_file = input_file
        chr_coverage_all =  utils.get_chr_coverage( fragments_file, [0,1000], chrlens, out=log )
        pickle.dump( chr_coverage_all, open( output_file[0],'w' ) , pickle.HIGHEST_PROTOCOL )
        
    return chr_coverage_all    

def insert_coverage_per_gene( input_file, output_file, annotation, 
                             mode='tss', blacklist = '', whitelist = '' ):
    
    
    '''
    Calculates insert coverage around TSS of each gene in the annotation. The interval (-1000,1000) around 
    the TSS is binned into 200 bp intervals and the depth normalized number of insertions is calculated for every bin. 
    Alternatively the annotation can contain centers of the peaks and then the function calculates the number of insertions 
    in -250..250 bin around the center of the peak. ('DHS' mode) 
    
    The depth normalization is done so that either the total number of insertions in the whitelist intervals is constant (if whitelist is given)
    or the total number of insertions everywhere outside the blacklist is constant. 
    
    Input: 
        input_file - .pickle file (output of get_coverage_vectors_atac_seq) with insertion densities    
        output_file - TSV file with rows corresponding to the names of the genes and columns to the bins around TSS
        annotation - gff or bed file with the annotation. **Every gene should contain a single interval** with the beginning 
                    corresponding to the TSS. 
        mode - TSS or DHS corresponds to if the annotation contains gene intervals or just single points of interest to calculate coverage around. 
        blacklist - TSV file with regions that should be removed from normalizing the total number of insertions 
        whitelist - TSV file with regions that should only be included when normalizing the number of insertions 
    
    Notes: 
        1. blacklist and whitelist parameters are mutually exclusive
        2. blacklist/whitelist lines are in the format: chromosome start end or just chromosome. Start and end points are included and 1-based. 
        
    Output: 
        None (output file is generated)
    '''
    
  
    dataframe = pd.DataFrame()
    if mode == 'tss':
        ranges = [ ( x-200,x ) for x in range( -800, 1200, 200 ) ]
    elif mode == 'dhs':
        ranges = [ ( -250,249 ) ]        
    
    for rg in ranges:
        
        if annotation.endswith('gff'):
            intervals = [ genomic_utils.get_interval_from_offset(x, rg, True) for x in 
                         genomic_utils.parse_gff_line_gen(open(annotation)) ]
            
        elif annotation.endswith('bed') or annotation.endswith('narrowPeak'): 
            intervals = [ genomic_utils.get_interval_from_offset(x, rg, True) for x in 
                         genomic_utils.parse_bed_line_gen(open(annotation)) ]
            
        icov = utils.get_ins_coverage(input_file, intervals, blacklist=blacklist,     
                                      whitelist=whitelist)
        for x in icov:
            icov[x] = icov[x].sum()
            
        dataframe[str(rg[0])] = pd.Series(icov)
        
    dataframe.to_csv(output_file, sep='\t')
    

def generate_tracks( input_file, output_files, chrlens, 
                    blacklist_file = '', whitelist_file = '', 
                    smoothing_profile='hanning', smoothing_window=1):
    '''
    Generates (potentially smoothed) depth normalized bigwig tracks
    Inputs:
        input_file - pickle file of insertion coverages (output of get_chr_coverage)
        output_files -  output file, logfile, errfile
        chrlen_file - TSV file with columns: chromosome, chromosome length
        
        smoothing_profile - profile of smoothing window (flat, hanning, hamming, bartlett, blackman)
        smoothing_window - size of smoothing window (default: 1, no smoothing)
        blacklist - TSV file with regions that should be removed from normalizing the total number of insertions 
        whitelist - TSV file with regions that should only be included when normalizing the number of insertions 
    
    Notes: 
        1. blacklist and whitelist parameters are mutually exclusive
        2. blacklist/whitelist lines are in the format: chromosome start end or just chromosome. Start and end points are included and 1-based. 

    
    Assumes: 
        bedGraphToBigWig (from UCSC userApps) command is on the $PATH
        
    '''
    
    input_data = pickle.load(open(input_file))
    logfile = open(output_files[1],'w',0)
    errfile = open(output_files[2],'w',0)
    logfile.write('Writing tracks\n')
    
        
    data=input_data
    if not blacklist_file :
        if not whitelist_file:
            norm_factor = utils.get_norm_factor(data)
        else: 
            whitelist = open( whitelist_file )
            norm_factor = utils.get_norm_factor( data, whitelist = whitelist)
            whitelist.close()
    else: 
        blacklist = open(blacklist_file)            
        norm_factor = utils.get_norm_factor(data, blacklist = blacklist)
        blacklist.close()
    for x in data:
        data[x] *= norm_factor
    logfile.write("Writing %s bedGraph\n"%output_files[0])
    f = tempfile.NamedTemporaryFile(dir= 
                                    os.environ.get('LOCAL_SSD','/tmp/'),
                                                  delete=False)
    
    chrvec = utils.get_chr_coverage_vector(data,
                            smoothing_window=smoothing_window,
                            smoothing_profile=smoothing_profile,
                            logfile=logfile)
    f.close()
    utils.write_chr_coverage_vectors(chrvec, f.name,mode='bedGraph')
    logfile.write ('Executing' + ' '.join(['bedGraphToBigWig', f.name, chrlens,
                                           output_files[0]]) + '\n')    
    try: 
        subprocess.check_call(['bedGraphToBigWig', f.name, chrlens, 
                               output_files[0]], stdout=logfile, stderr=errfile)
    except OSError:  
        print >>sys.stderr, 'bedGraphToBigWig not found'
        sys.exit(-1)
    except subprocess.CalledProcessError:
        print >>sys.stderr, ' '.join(['bedGraphToBigWig', f.name, chrlens,
                                      output_files[0]]), 'crashed'
        sys.exit(-1)
    
    f.unlink(f.name)
    logfile.close()
    errfile.close()

