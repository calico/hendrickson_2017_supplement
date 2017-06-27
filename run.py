#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:39:43 2017

@author: ilyasoifer
"""

from varutils import mkdir_p, mkname
import argparse, sys
from os.path import isfile


import pipeline

def is_valid_file(parser, arg):
    if arg and not exists(arg) :
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle



def parse_args(  ):
    parser = argparse.ArgumentParser(description=
                                        'Example of ATACseq processing', 
                                        prog='atacseq_processing')
    parser.add_argument("input_file", help="Unsorted BAM file", 
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("output_dir", help="Destination directory", type=str)
    parser.add_argument("chrlens_file", help="Chromosome length file",  
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("annotation_file", help="Annotation file", 
                        type=lambda x: is_valid_file(parser, x))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--whitelist", help="Name of whitelist file", default='', 
                       type=lambda x: is_valid_file(parser, x))
    group.add_argument("--blacklist", help="Name of blacklist file", default='',
                       type=lambda x: is_valid_file(parser, x))
    return parser.parse_args()


params = parse_args()
input_file = params.input_file
output_dir = params.output_dir
chrlens_file = params.chrlens_file
mkdir_p( output_dir )



########## STEP1: BAM to BED ################
intervals_file = mkname( output_dir, input_file, '.bam','.bed')
log_file = mkname( output_dir, input_file, '.bam','.bed.log')
err_file = mkname( output_dir, input_file, '.bam','.bed.err')


print >>sys.stdout, "Converting BAM (%s) to BED (%s)"%( input_file, intervals_file )
sys.stdout.flush()
pipeline.bam2intervals( input_file, [ intervals_file, log_file, err_file ], chrlens_file )

########## STEP2: BED to CVG vectors ########
cvg_file = mkname( output_dir, input_file, '.bam','.cvg.pickle')
log_file = mkname( output_dir, input_file, '.bam','.cvg.log')


print >>sys.stdout, "Converting BED (%s) to coverages (%s)"%( intervals_file, cvg_file )
pipeline.get_coverage_vectors_atac_seq( intervals_file, [ cvg_file, log_file ], chrlens_file )


########## STEP3: CVG of annotation #########
annotation_file = params.annotation_file
whitelist = params.whitelist
blacklist = params.blacklist
annotation_cvg = mkname( output_dir, input_file, '.bam','.annotation.tsv' )

print >>sys.stdout, "Converting CVG (%s) to annotation coverages (%s)"%( cvg_file, annotation_file )

pipeline.insert_coverage_per_gene(cvg_file, annotation_cvg, annotation_file, whitelist=whitelist,
                                  blacklist = blacklist, mode='tss')


########## STEP4: CVG to bigWig track #######
bigwig_file = mkname( output_dir, input_file, '.bam','.bigWig')
log_file = mkname( output_dir, input_file, '.bam','.bigWig.log')
err_file = mkname( output_dir, input_file, '.bam','.bigWig.err')

print >>sys.stdout, "Converting CVG (%s) to bigWig track (%s)"%( cvg_file, bigwig_file )

pipeline.generate_tracks(cvg_file, [ bigwig_file, log_file, err_file ], chrlens_file, 
                         blacklist_file = blacklist, whitelist_file=whitelist, 
                         smoothing_profile='hanning', smoothing_window=1)

