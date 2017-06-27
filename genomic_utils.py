#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:18:42 2017

@author: Ilya
"""

import numpy as np
import types, re

class Position:
    def __init__( self, tup, weight=None ):
        '''
        tup is either tab separated chr\tpos\t\tlength\t\tstrand
        or tuple (chr, pos, length, strand)
        weight could be string ( converted to float then divided by 100), int (then divided by 100) or float
        '''
        if type( tup ) == type(''):
            tup = tup.strip().split('\t')
            self.chromosome = tup[0]
            self.pos        = int(tup[1])
            self.length     = int(tup[3])
            self.strand     = tup[5]
            self.gene = hash(self)
        else:
            if len(tup) == 4:
                self.chromosome, self.pos, self.length, self.strand = tup
                self.gene = hash(self)
            elif len(tup) == 5 : 
                self.chromosome, self.pos, self.length, self.strand, self.gene = tup
            if weight:
                if type(weight) == types.StringType or type(weight) == types.IntType: 
                    self.weight = float(weight)/100
                elif type(weight) == types.FloatType:
                    self.weight = weight
                else :
                    raise(RuntimeError("weight should be numeric or convertible to numeric"))
            
    
    def __str__(self)  :
        return 'Position: %s(%s):%d'%(self.chromosome, self.strand, self.pos)

    def __hash__(self):    
        return hash( ( self.chromosome, self.pos, self.length, self.strand ))
        
        
class Interval:
    def __init__ ( self, tup, ignore_strand_flag = None ):
        '''
        tup can be 
        1. String, then it is split to >= 6 fields: (chromosome, start, end, .., .., strand)
        2. Position
        3. Interval (that may have Gene interval)
        4. Tuple (chromsome, start, end, strand, [gene])
        The interval is always zero-based right open. 
        '''
        if type(tup) == type(''):
            tup = tup.strip().split('\t')
            self.gene       = None
            self.chromosome = tup[0]
            self.start      = int( tup[1] )
            self.end        = int( tup[2] )
            self.strand     = int( tup[5] )
        elif isinstance( tup, Position):
            self.chromosome = tup.chromosome
            self.start = tup.pos - tup.length/2
            self.end   = tup.pos + tup.length/2 + 1
            self.strand = tup.strand
        elif isinstance( tup, Interval ):
            self.chromosome = tup.chromosome
            self.start      = tup.start
            self.end        = tup.end
            self.strand     = tup.strand
            self.ignore_strand = tup.ignore_strand
            if hasattr(tup,'gene'):
                self.gene       = tup.gene
        else:
            self.chromosome = tup[0]
            self.start      = tup[1]
            self.end        = tup[2]
            self.strand     = tup[3]
            if len(tup)>=5 :
                self.gene       = tup[4]
        if type(ignore_strand_flag) == type(False) :
            self.ignore_strand = ignore_strand_flag
        if not ignore_strand_flag :
            self.ignore_strand = False
        if not hasattr(self,'ignore_strand'):
            raise RuntimeError('Interval must have ignore_strand flag')
            
        
    def __len__( self ):
        if not hasattr(self, 'end'):
            return np.NaN
        return self.end - self.start 
    
    def __hash__( self ):
        return hash( ( self.chromosome, self.start, self.end, self.strand ))

    def get_pos( self, x ):
        if x < len(self) and x >= 0:
            if self.ignore_strand or self.strand == '+':
                if hasattr(self,'gene'):
                    return Position((self.chromosome, self.start + x , 1, self.strand, self.gene))
                else:
                    return Position((self.chromosome, self.start + x, 1, self.strand))
            else:
                if hasattr(self,'gene'):
                    return Position((self.chromosome, self.end - 1 - x , 1, self.strand, self.gene))
                else:
                    return Position((self.chromosome, self.end - 1 - x, 1, self.strand))
        else:
            raise IndexError("Position out of range")
        
    def to_generator(self ):
        '''
        generates a tuples of consecutive positions in the interval
        '''
        chromosome = self.chromosome
        start      = self.start
        end        = self.end
        strand     = self.strand
        if strand  == '+':
            pos_place = -1
        else : 
            pos_place = end - start 
        for i in range(start, end):
            if strand == '+' :
                pos_place +=1 
            else: 
                pos_place -=1 
            yield( chromosome, i, strand, pos_place  )
    def compare_to_pos( self , pos ):
        if pos.chromosome < self.chromosome:
           
            return -1
        if pos.chromosome > self.chromosome:
            #print "there"
            return 1
        #print pos.chromosome, self.chromosome
        if pos.pos < self.start :
            return -1
        if pos.pos >= self.end:
            
            return 1
        if not self.ignore_strand:
            
            if pos.strand < self.strand:
                return -1
            if pos.strand > self.strand: 
                return 1
        
        return 0
    
    def __cmp__(self, other):
        return ( self.chromosome, self.start, self.strand) < ( other.chromosome, other.start, other.strand )
        
    def __str__(self):
        return '%s(%s) %d-%d'%(self.chromosome, self.strand, self.start, self.end)
        
    def pos_in_interval( self, position ):
        if not self.ignore_strand:
            if self.strand != position.strand:
                return -1
        if self.chromosome != position.chromosome:
            return -1
        if self.start > position.pos or self.end <= position.pos:
            return -1
        if self.strand == '+':
            return position.pos - self.start
        else: 
            return self.end - 1 - position.pos #CAREFUL: check how intervals are generated

    def dist_to_interval( self, position ):
        if self.pos_in_interval( position ) > 0:
            return 0
        elif self.chromosome != position.chromosome:
            return np.Inf
        else:
            to_start = self.start - position.pos
            to_end   = position.pos - self.end-1
            if self.strand == '+':
                if abs(to_start) <= abs(to_end) :
                    return -to_start
                else: 
                    return to_end
            else: 
                if abs(to_start) <= abs(to_end) :
                    return to_start
                else:
                    return -to_end
            
    def to_bed( self, score = None, rgb = None ):
        '''
        Generates a bed representation of the interval (optionally with score and rgb)
        '''
        strg = '%s\t%d\t%d'%(self.chromosome, self.start, self.end)
        if hasattr(self,'gene'):
            strg += "\t%s"%self.gene
        else:
            strg += "\t%d"%hash(self)
        if score:
            strg += "\t" + str(score)
        else:
            strg += "\t"
        if self.ignore_strand:
            strg += "\t"
        else:
            strg += "\t" + self.strand
        
        strg +="\t\t\t"
        if rgb:
            strg += "%d,%d,%d"%rgb
        return strg
        
    def dist( self, other ):
        '''
        Minimal distance between two intervals
        0 if there is an intersection
        '''
        if self.chromosome != other.chromosome:
            return np.inf
        if not self.ignore_strand or not other.ignore_strand and self.strand != other.strand :
            return np.inf
        spos = other.get_pos(0)
        epos = other.get_pos( len( other ) - 1 )
        res = min( self.dist_to_interval( spos ), self.dist_to_interval( epos ) )
        spos = self.get_pos(0)
        epos = self.get_pos( len( self ) - 1 )
        res = min(res, other.dist_to_interval( spos ) )
        res = min(res, other.dist_to_interval( epos ) )
        return res
        
get_chr_lens = lambda chrlen_file: dict([ (x.split()[0],int(x.split()[1])) for x in open(chrlen_file) ] )

def parse_gff_line( strg, feature_type_filter = 'CDS', name_regexp=r'.*Name=([a-zA-Z0-9-_]*);' ):
    fields = strg.split('\t')
    #print len(fields)
    if len(fields) < 7: 
        return
    chrom        = fields[0]
    feature_type = fields[2]
    if feature_type_filter and feature_type != feature_type_filter:
        return
    start        = int(fields[3]) - 1 #GFF is one-based
    end          = int(fields[4]) # GFF is both sides closed
    strand       = fields[6] 
    name         = re.match(name_regexp, fields[8])
    if name:
        name = name.groups(0)[0]
    else:
        name = ''
    return ( chrom, start, end, strand, name, 1 )


def parse_gff_line_gen(gff_iterator, feature_type_filter = 'CDS', name_regexp=r'.*Name=([a-zA-Z0-9-_]*);'):
    
    for line in gff_iterator:
        if not line or line.startswith('#'):
            continue
        tup = parse_gff_line( line, feature_type_filter, name_regexp )
        if not tup:
            continue
        yield tup 
    return

def parse_bed_line_gen(bed_iterator):
    '''
    returns ( chrom, start, end, strand, name, score ) for each line
    '''
    
    for line in bed_iterator:
        if not line or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom        = fields[0]
        start        = int(fields[1]) #BED is zero based
        end          = int(fields[2]) # BED is half-open
        if len( fields ) < 6 :
            strand = '+'
        else:
            strand       = fields[5] 
        if len(fields) < 4:
            name = 'Noname'
        else:
            name         = fields[3]
        if len( fields ) >= 5 and fields[4]:
            score = float(fields[4])/100
        else:
            score = 1.0
        yield ( chrom, start, end, strand, name, score )
    return

def get_interval_from_offset(interval, offset, ignore_strand):
    interval     = Interval(interval, ignore_strand)
    new_interval = Interval(interval)
    
    if new_interval.strand == '+':
        new_interval.start = interval.start + offset[0]
        new_interval.end   = interval.start + offset[1]
    else:
        new_interval.start = interval.end  - offset[1]
        new_interval.end   = interval.end  - offset[0]
    return new_interval
