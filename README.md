# Supplementary code for Henderson et al., 2017 "Simultaneous profiling of DNA accessibility and gene expression dynamics with ATAC-seq and RNA-seq" 

## Summary 

This is the code for reproducing the bioinformatic analyses of the yeast ATACseq data (section 3.3). We assume that the library was demultiplexed and aligned to the genome as described in sections (3.3.1-3.3.2). Analysis of RNAseq data is described in detail elsewhere (reference to the nature protocols of cufflinks). 

The processing is done by `running atacseq_processing.py` as described below. The input is the aligned __unsorted__ BAM file and the output is the bigWig track of normalized insertion densities, python `pickle` file that contains the insertion coverage for every nucleotide in the genome and tab separated file containgin the coverage of insertions in 200 bp bins spanning -1000..1000 bases around the transcription start sites. 

## Requirements 

* python2.7
* numpy
* scipy
* pandas
* pysam (>= 0.18)
* UCSC command line bioinformatic utils (a.k.a. kentUtils)

## Installation

### Clone git repository

    mkdir atacseq_analysis
    git clone ??? atacseq_analysis

### Download example data
    cd atacseq_analysis
    mkdir data
    cd data 
    download from bucket

## Usage

    usage: atacseq_processing.py [-h]
                                 [--whitelist WHITELIST | --blacklist BLACKLIST]
                                 input_file output_dir chrlens_file
                                 annotation_file

    Example of ATACseq processing

    positional arguments:
      input_file            Unsorted BAM file
      output_dir            Destination directory
      chrlens_file          Chromosome length file
      annotation_file       Annotation file

    optional arguments:
      -h, --help            show this help message and exit
      --whitelist WHITELIST
                            Name of whitelist file
      --blacklist BLACKLIST
                            Name of blacklist file
                            
## Example 

    python atacseq_processing.py gcn processed_data data/SacCer3.fasta.len data/annotation.bed 
                            --blacklist data/blacklist.tsv 
