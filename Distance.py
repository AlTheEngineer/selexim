#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import sys, time
import random
import linecache
from itertools import izip, imap, islice
import operator
from collections import OrderedDict
import numpy as np
from scipy import stats
import utils
from utils import apt_loopFinder
#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance

class Distance:
#This function takes the sequences of two loop regions and returns
#their Lavenshtein distance
#Input: str(), str()
#Output: int()
    def lavenshtein_func(self , loop1, loop2):
        #if statement for switching loops
        if len(loop1) < len(loop2):
            return self.lavenshtein_func(loop2, loop1)
        if len(loop2) == 0:
            return len(loop1)
        #separate individual nucleotides in each loop
        loop1 = np.array(tuple(loop1))
        loop2 = np.array(tuple(loop2))
        prev_row = np.arange(loop2.size + 1)
        for nt in loop1:
            curr_row = prev_row +1
            curr_row[1:] = np.minimum(curr_row[1:], np.add(prev_row[:-1], loop2 != nt))
            curr_row[1:] = np.minimum(curr_row[1:], curr_row[0:-1] + 1)
            prev_row = curr_row
        loop2_dist = prev_row[-1]
        return loop2_dist


#This function takes two sequences of equal length and returns
#their Hamming distance
#Input: str(), str()
#Output: int()
    def hamming_func(self, seq1, seq2):
        len(seq1) == len(seq2)
        ne = operator.ne
        seq2_dist = sum(imap(ne, seq1, seq2))
        return seq2_dist

#This function takes the secondary structure of the reference aptamer
#and an arbitrary sequence and returns
#their Base-pair distance
#Input: str(), str()
#Output: int()
    def bp_func(self, seq1_struct, seq2):
        seq2_struct = fold(seq2)[0]
        seq2_dist = bp_distance(seq1_struct, seq2_struct)
        return seq2_dist

#This function takes the sequence, loop region and secondary structure of the reference aptamer
#and an arbitrary sequence and their lengths and returns the Loop-based distance
#Input: str(), str(), str(), str(), int()
#Output: int()
    def loop_func(self, seq1, seq1_struct, seq1_loop, seq2, seqLength):
        #compute secondary structure of sequence
        seq2_struct = fold(seq2)[0]
        base = None
        baseIdx = 0
        #find a 3' paired nucleotide
        while(base != ')' and baseIdx < seqLength-1):
            base = seq2_struct[baseIdx]
            baseIdx += 1
        if(baseIdx == seqLength-1):
            while(base != '(' and baseIdx > 0):
                base = seq2_struct[baseIdx-1]
                baseIdx -= 1
            if(baseIdx == 0):
                #sequence doesnt have a loop
                seq2_loop = seq2
            else:
                #sequence loop is dangling end
                seq2_loop = seq2[baseIdx:]
        else:
            #sequence has a loop
            loop_end = baseIdx-1
            while(base != '('):
                baseIdx -= 1
                base = seq2_struct[baseIdx-1]
            #grab loop
            seq2_loop = seq2[baseIdx:loop_end]
        #compute Lavenshtein distance
        seq2_loopDist = self.lavenshtein_func(seq1_loop, seq2_loop)
        #compute BP distance
        seq2_bpDist = bp_distance(seq1_struct, seq2_struct)
        #sum distances
        seq2_dist = int(seq2_loopDist + seq2_bpDist)
        return seq2_dist
#TEST
#d = Distance()
#seq_dist = d.loop_func(he4_seq, he4_struct, ex_seq)


#This function takes the sequence, loop region and secondary structure of the reference aptamer
#and an arbitrary sequence and their lengths and returns the component Lavenshtein and BP
#distances
#Input: str(), str(), str(), str(), int()
#Output: int(), int()
    def loop_components_func(self, seq1, seq1_struct, seq1_loop, seq2, seqLength):
        seq2_struct = fold(seq2)[0]
        base = None
        baseIdx = 0
        while(base != ')' and baseIdx < seqLength-1):
            base = seq2_struct[baseIdx]
            baseIdx += 1
        if(baseIdx == seqLength-1):
            while(base != '(' and baseIdx > 0):
                base = seq2_struct[baseIdx-1]
                baseIdx -= 1
            if(baseIdx == 0):
                seq2_loop = seq2
            else:
                seq2_loop = seq2[baseIdx:]
        else:
            loop_end = baseIdx-1
            while(base != '('):
                baseIdx -= 1
                base = seq2_struct[baseIdx-1]
            seq2_loop = seq2[baseIdx:loop_end]
        seq2_loopDist = self.lavenshtein_func(seq1_loop, seq2_loop)
        seq2_bpDist = bp_distance(seq1_struct, seq2_struct)
        return seq2_loopDist, seq2_bpDist

#This function takes a sequence and its length and computes its bias score
#Input: str(), int()
#Output: float()
    def bias_func(self, seq, seqLen):
        pyrNum = 0
        for nt in seq[:-1]:
            if(nt == 'C') or (nt == 'T'):
                pyrNum += 1 #increment no. of pyrimidines
        biasScore = 0.1*(2*pyrNum - seqLen)/seqLen #compute bias
        return biasScore
