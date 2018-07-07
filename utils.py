#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import sys, time
import random
import linecache
import math
from math import factorial
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np
from scipy.misc import comb
from scipy.interpolate import spline
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import Distance

#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance, svg_rna_plot



#This computes the number of total and unique sequences
#present in the pool
def seqNumberCounter(seqPool):
    totalSeqNum = int(0)
    uniqSeqNum = int(0)
    for seqIdx in seqPool:
        totalSeqNum += seqPool[seqIdx][0]
        uniqSeqNum += 1
    return int(totalSeqNum), int(uniqSeqNum)
#This computes the binomial coefficient (not used)
def binomCoeff(n, k):
    binom = factorial(n)/(factorial(k)*factorial(n-k))
    return binom

#This converts an array of probabilities into a 
#discrete probability distribution
def convert_to_distribution(x, y, distName):
    xDist = stats.rv_discrete(name=distName, values=(x, y))
    return xDist

#This finds the loop region in a given sequence
def apt_loopFinder(apt_seq, apt_struct, seqLength):
    base = None
    baseIdx = 0
    while(base != ')' and baseIdx < seqLength):
        base = apt_struct[baseIdx]
        baseIdx += 1
    if(baseIdx == seqLength):
        while(base != '('and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:]
            return apt_loop
    else:
        loop_end = baseIdx-1
        while(base != '(' and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq[:loop_end]
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:loop_end]
            return apt_loop

# Add method for computing the L1 norm 

#This constructs a discrete probability distribution to be used for 
#sampling during the selection step
def rvd(X, X_sum, distName):
    seqIdxs = np.zeros(X.shape[0])
    probs = np.zeros(X.shape[0])
    for i, seq in enumerate(X):
        seqIdxs[i] = i
        probs[i] = seq[1]/X_sum
    dist = stats.rv_discrete(name=distName, values=(seqIdxs, probs))
    return dist
#This calculates the overall average bias for total and unique sequences
#for the given round
def bias_avg(seqFile, seqLength):
    bias = 0
    w_bias = 0
    totalSeqs = 0
    uniqSeqs = 0
    with open(seqFile, 'r') as f:
        for line in f:
            row = line.split()
            seq = row[0]
            bias += d.bias_func(seq, seqLength)
            w_bias += int(row[1])*d.bias_func(seq, seqLength)
            totalSeqs += int(row[1])
            uniqSeqs += 1
    avg_bias = bias/uniqSeqs
    weighted_avg_bias = w_bias/totalSeqs
    return weighted_avg_bias, avg_bias

#This calculates the average bias for each affinity group in the given round
def bias_avg_per_dist(seqFile, seqLength):
    bias_per_dist = np.zeros(seqLength+5)
    w_bias_per_dist = np.zeros(seqLength+5)
    totalSeqs_per_dist = np.zeros(seqLength+5)
    uniqSeqs_per_dist = np.zeros(seqLength+5)
    with open(seqFile, 'r') as f:
        for line in f:
            row = line.split()
            #grab sequence
            seq = row[0]
            #grab frequency
            freq = int(row[1])
            #grab distance
            dist = int(row[2])
            uniqSeqs_per_dist[dist] += 1
            totalSeqs_per_dist[dist] += freq
            #add bias to corresponding dist
            bias_per_dist[dist] += d.bias_func(seq, seqLength)
            w_bias_per_dist[dist] += freq*d.bias_func(seq, seqLength)
    #calculate averages
    for dist in xrange(seqLength+5):
        if(uniqSeqs_per_dist[dist] > 0):
            bias_per_dist[dist] /= uniqSeqs_per_dist[dist]
            w_bias_per_dist[dist] /= totalSeqs_per_dist[dist]
    return w_bias_per_dist[:seqLength+1], bias_per_dist[:seqLength+1]


#bias_avg("window_R14", 20)

#This plots the theoretical distribution of Hamming distances
def seq_div_hamm(seqLength, alphabetSet):
    uniqSeqNum_per_dist = np.zeros(seqLength+1)
    for h in xrange(seqLength+1):
        uniqSeqNum_per_dist[h] = (len(alphabetSet)-1)**(h)*comb(seqLength, h)
    hammDistAxis = np.linspace(0, seqLength, seqLength+1)
    hammDistAxis_smooth = np.linspace(0, seqLength, 200)
    uniqSeqNum_smooth = spline(hammDistAxis, uniqSeqNum_per_dist, hammDistAxis_smooth)
    fig, ax = plt.subplots(1,1)
    ax.plot(hammDistAxis_smooth, uniqSeqNum_smooth)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fig.text(0.5, 0.95, 'Unique Sequences', ha='center')
    fig.text(0.5, 0.04, 'Hamming Distance', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    fig.savefig("SELEX_Analytics_seqDiv_20nt", format='pdf')
    return 0
#seq_div_hamm(20, 'ACGT')

#This plots the range of distances for the Hamming, BP, and
#Loop-based metric. The scale defines the number of samples
#to use to construct the distribution 
def distance_range(scale, ref_seq, seqLength, alphabetSet):
    ref_struct = fold(ref_seq)[0]
    ref_loop = apt_loopFinder(ref_seq, ref_struct)
    hamm_dist_array = np.zeros(int(seqLength*1.5))
    bp_dist_array = np.zeros(int(seqLength*1.5))
    loop_dist_array = np.zeros(int(seqLength*1.5))
    randIdxs = random.randint(0, 4**(20)-1, size=scale)
    for i in xrange(scale):
        randIdx = randIdxs[i]
        randSeq = apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
        randHammDist = d.hamming_func(randSeq, ref_seq)
        randbpDist = d.bp_func(ref_struct, randSeq)
        randLoopDist = d.loop_func(ref_seq, ref_struct, ref_loop, randSeq, seqLength)
        hamm_dist_array[randHammDist] += 1
        bp_dist_array[randbpDist] += 1
        loop_dist_array[randLoopDist] += 1
    for dist in xrange(int(seqLength*1.5)):
        hamm_dist_array[dist] /= scale
        bp_dist_array[dist] /= scale
        loop_dist_array[dist] /= scale
    fig, axis = plt.subplots(1,1)
    distAxis = np.linspace(0, int(seqLength+9), int(seqLength+10))
    distAxis_smooth = np.linspace(0, int(seqLength+9), 200)
    hamm_dist_smooth = spline(distAxis, hamm_dist_array, distAxis_smooth)
    bp_dist_smooth = spline(distAxis, bp_dist_array, distAxis_smooth)
    loop_dist_smooth = spline(distAxis, loop_dist_array, distAxis_smooth)
    axis.plot(distAxis_smooth, hamm_dist_smooth, label='Hamming')
    axis.plot(distAxis_smooth, bp_dist_smooth, label='Base-Pair')
    axis.plot(distAxis_smooth, loop_dist_smooth, label='Loop')
    axis.set_xlim([0, 25])
    axis.set_ylim([0, 0.4])
    axis.legend()
    fig.text(0.5, 0.04, 'Distance', ha='center')
    fig.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
    fig.text(0.5, 0.95, 'Distance Distributions', ha='center')
    fig.savefig("SELEX_Analytics_distance_distributions", format='pdf')
    return hamm_dist_array

#Identify the aptamer candidate with the highest frequency in the final pool
#outputs the candidate's structure in .svg format
def aptamer_structs(fileNames, seqLength, roundNum, rounds='final'):
    if(rounds == 'final'):
        top_seq_info = [0,0]
        with open(fileNames+"_R"+str(roundNum), 'r') as f:
            for line in f:
                row = line.split()
                seq = str(row[0])
                count = int(row[1])
                if(count > top_seq_info[1]):
                    top_seq_info[0] = seq
                    top_seq_info[1] = count
        with open(fileNames+"_R"+str(roundNum)+"_topstructure_info", 'w') as f:
            seq = top_seq_info[0]
            seq_struct = fold(seq)[0]
            seq_mfe = fold(seq)[1]
            seq_count = top_seq_info[1]
            f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\n')
        svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(roundNum)+"_topstructure.svg")
        return 0
    elif(rounds == 'all'):
        top_seqs_info = []
        for rnd in xrange(roundNum):
            with open(fileNames+"_R"+str(rnd+1), 'r') as f:
                for line in f:
                    row = line.split()
                    seq = str(row[0])
                    count = int(row[1])
                    if(count > top_seq_info[1]):
                        top_seqs_info.append([seq,count])
        with open(fileNames+"_R"+str(roundNum)+"_topstructures_info", 'w') as f:
            for rnd in xrange(roundNum):
                seq = top_seqs_info[rnd][0]
                seq_struct = fold(seq)[0]
                seq_mfe = fold(seq)[1]
                seq_count = top_seqs_info[rnd][1]
                f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\n')
                svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(rnd+1)+"_topstructure.svg")
        return 0
    else:
        print("invalid option for string varible rounds. Exiting...")

#aptamer_structs('he4_hamm_small', 20, 40, 'final')


#identify the aptamer candidates with the highest affinity in the final pool
#outputs the structure of the candidate in .svg format
def aptamer_structs_aff(fileNames, seqLength, roundNum, rounds='final'):
    if(rounds == 'final'):
        top_seq_info = [0,0,np.infty]
        with open(fileNames+"_R"+str(roundNum), 'r') as f:
            for line in f:
                row = line.split()
                seq = str(row[0])
                count = int(row[1])
                dist = int(row[2])
                if(dist < top_seq_info[2]):
                    top_seq_info[0] = seq
                    top_seq_info[1] = count
                    top_seq_info[2] = dist
        with open(fileNames+"_R"+str(roundNum)+"_affstructure_info", 'w') as f:
            seq = top_seq_info[0]
            seq_struct = fold(seq)[0]
            seq_mfe = fold(seq)[1]
            seq_count = top_seq_info[1]
            seq_dist = top_seq_info[2]
            f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\t'+str(seq_dist)+'\n')
        svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(roundNum)+"_affstructure.svg")
        return 0
    elif(rounds == 'all'):
        top_seqs_info = []
        for rnd in xrange(roundNum):
            with open(fileNames+"_R"+str(rnd+1), 'r') as f:
                for line in f:
                    row = line.split()
                    seq = str(row[0])
                    count = int(row[1])
                    dist = int(row[2])
                    if(dist > top_seq_info[2]):
                        top_seqs_info.append([seq,count,dist])
        with open(fileNames+"_R"+str(roundNum)+"_affstructures_info", 'w') as f:
            for rnd in xrange(roundNum):
                seq = top_seqs_info[rnd][0]
                seq_struct = fold(seq)[0]
                seq_mfe = fold(seq)[1]
                seq_count = top_seqs_info[rnd][1]
                seq_dist = top_seqs_info[rnd][2]
                f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)++'\t'+str(seq_dist)+'\n')
                svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(rnd+1)+"_affstructure.svg")
        return 0
    else:
        print("invalid option for string varible rounds. Exiting...")

#aptamer_structs_aff('he4_bp_small', 20, 40, 'final')


