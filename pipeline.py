#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:27:58 2017

@author: Ilya
"""
import utils, genomic_utils
from copy import deepcopy
import pickle, tempfile, os, subprocess, sys
import numpy as np
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
        input_file - bed file of inserts (output of bam2intervals)
        output_file[0] - pickle file with dictionary (per-chromosome) of insert coverages
        output_file[1] - log file
        chrlen_file - TSV file with columns: chromosome, chromosome length
    
    Output
    See also: 
        utils.get_chr_coverage
        utils.get_chr_coverage_vector
    '''
    with open( output_file[1], 'w' ) as log :
        fragments_file = input_file
        chr_coverage_all =  utils.get_chr_coverage( fragments_file, [0,1000], chrlens, out=log )
        pickle.dump( chr_coverage_all, open( output_file,'w' ) , pickle.HIGHEST_PROTOCOL )
        
    return chr_coverage_all    

def insert_coverage_per_gene( input_file, output_file, annotation, part = 'all', 
                             mode='tss', blacklist = '', whitelist = '' ):
    '''
    Calculates insert coverage per gene in the annotation
    input_file - [.pickle]
    output_file - []
    
    annotation - gff or bed 
    part - should be high or low or all (default: all)
    '''
    data = pickle.load(open(input_file[0]))[part]
    dataframe = pd.DataFrame()
    if mode == 'tss':
        ranges = [ (x-200,x) for x in range( -800, 1200, 200 ) ]
    elif mode == 'dhs':
        ranges = [ (-250,249) ]        
    
    for rg in ranges:
        if annotation.endswith('gff'):
            intervals = [ genomic_utils.get_interval_from_offset(x, rg, True) for x in 
                         genomic_utils.parse_gff_line_gen(open(annotation)) ]
        elif annotation.endswith('bed') or annotation.endswith('narrowPeak'): 
            intervals = [ genomic_utils.get_interval_from_offset(x, rg, True) for x in 
                         genomic_utils.parse_bed_line_gen(open(annotation)) ]
        icov = utils.get_ins_coverage(data, intervals, blacklist=blacklist, 
                                               normlist=whitelist)
        for x in icov:
            icov[x] = icov[x].sum()
        dataframe[str(rg[0])] = pd.Series(icov)
    dataframe.to_csv(output_file[0], sep='\t')

def generate_tracks( input_file, output_files, chrlens, smoothing_window=1, 
                    blacklist_file = '', whitelist_file = '', 
                    smoothing_profile='hanning'):
    '''
    Generates (potentially smoothed) depth normalized bigwig tracks
    Inputs:
    input_file - list of single pickle file
    output_files: output for low, high and all fragments, logfile, errfile
    chrlens file 
    size of smoothing window (default 1) 
    '''
    input_file = input_file[0]
    insertions_parts = ['low','high','all']
    input_data = pickle.load(open(input_file))
    logfile = open(output_files[3],'w',0)
    errfile = open(output_files[4],'w',0)
    logfile.write('Writing tracks\n')
    
    for p in range(3):
        if not input_data.has_key(insertions_parts[p]) or not output_files[p]:
            continue
        
        data=input_data[insertions_parts[p]]
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
        logfile.write("Writing %s bedGraph\n"%insertions_parts[p])
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
                                               output_files[p]]) + '\n')    
        try: 
            subprocess.check_call(['bedGraphToBigWig', f.name, chrlens, 
                                   output_files[p]], stdout=logfile, stderr=errfile)
        except OSError:  
            print >>sys.stderr, 'bedGraphToBigWig not found'
            sys.exit(-1)
        except subprocess.CalledProcessError:
            print >>sys.stderr, ' '.join(['bedGraphToBigWig', f.name, chrlens,
                                          output_files[p]]), 'crashed'
            sys.exit(-1)
        
        f.unlink(f.name)
    logfile.close()
    errfile.close()


def read_and_normalize_atac_seq( directory, filenames, labels ):
    '''
    Reads the estimtated coverages of atac_seq from time series/mutants
    done by seq_pipeline_utils.insert_coverage_per_gene, unites them and normalizes by the expression of the first labek
    directory - where the .annot.exp files are located 
    filenames - list of the files in order of labels
    returns
    datasets - each timepoint/mutant separately as pandas Dataframe
    results  - combined matrix ( as pandas dataframe, normalized values will be prefixed with Change)
    
    the numbers annot.exp should be depth normalized
    '''
    fullfiles = [ os.path.join(directory,x) for x in filenames ]
    datasets = {}
    for i in range( len( labels ) ):
        label = labels[ i ]
        datasets[ label ] = pd.read_csv( fullfiles[ i ], delim_whitespace = True, index_col = 0 )

    result = datasets[labels[0]]
    n_columns = len(result.columns)
    orig_names = list(result.columns)
    for label in labels[1:]: 
        result = result.join(datasets[label], rsuffix = '_' + label , how='outer')
    cols = list(result.columns)
    for i in range(n_columns):
        cols[i] = cols[i] + '_' + labels[0] 
    result.columns = cols
    for name in orig_names:
        for label in labels[1:]:
            result['_'.join(['Change',label,name][::-1])] = np.log2(result['_'.join([name,label])]) - np.log2(result['_'.join([ name,labels[0]])])
    cols = result.columns
    new_cols = [ '_'.join(x.split('_')[::-1]) for x in cols]
    result.columns = new_cols
    submat = np.array(result[ [ x for x in result.columns if 'Change' in x ]])
    submat = varutils.complete_nans( submat )
    result[ [  x for x in result.columns  if 'Change' in x ] ] = submat
    return datasets,result

