# Supplementary code for Hendrickson et al., "Simultaneous profiling of DNA accessibility and gene expression dynamics with ATAC-seq and RNA-seq" 

## Summary 

This is the code for reproducing the bioinformatic analyses of the yeast ATACseq data (section 3.3). We assume that the library was demultiplexed and aligned to the genome as described in sections 3.3.1-3.3.2. Analysis of RNAseq data is described in detail [elsewhere](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html). 

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
    git clone git@github.com:calico/hendrickson_2017_supplement.git atacseq_analysis

### Download example data

    cd atacseq_analysis
    wget https://hendrickson_2017_supplement.storage.googleapis.com/data.tar.gz
    tar xvfz data.tar.gz
    

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
                            
## Parameter explanation
   * input_file - unsorted BAM file with aligned reads
   * output_dir - directory to write the output (will be created if missing)
   * chrlens_file - TSV of format: chromosome   length
   * annotation_file - BED or GTF file with a single interval per gene. 
   * whitelist - if given, the insertion coverage will be normalized so that the total
    number of insertions in the whitelist is the same between samples.
   * blacklist - if given, the insertion coverage will be normalized so that the total number of insertions in the genome
   outside the blacklist is constant between samples. 
   
   Both whitelist and blacklist are TSV files with lines of the format `chromosome  start end`, where the intervals 1-based 
   and closed. Alternatively, also lines with only the chromosome name are also permitted. 
   
## Example 
    cd atacseq_analysis
    python atacseq_processing.py data/GCN4_t0_S1_L001_R1_001.mrg.bam processed_data \
                            data/SacCer3.fasta.len data/annotation.bed \
                            --blacklist data/blacklist.tsv 
