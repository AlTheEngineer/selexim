#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import time
import gc
from scipy import stats
from scipy.stats import norm
import numpy as np
from numpy.random import binomial as binom
from numpy.random import poisson
import matplotlib
#matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
from sklearn.mixture import GMM
from sklearn.mixture import VBGMM
from sklearn.mixture import DPGMM
from matplotlib import rcParams
from Mutation import Mutation

#Allow LaTeX in plots
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
#Turn off interactive mode for plots
plt.ioff()


#Initiate class
class Amplification(object):
    # add __init__ contructor here

    def randomPCR_with_ErrorsAndBias_FASTv2(self, slctdSeqs, 
                                     seqLength, pcrCycleNum, 
                                     pcrYld, errorRate, 
                                     aptamerSeqs, alphabetSet, distance):
        # initialize Mutation object from class
        mut = Mutation(seqLength=seqLength, errorRate=errorRate, 
                        pcrCycleNum=pcrCycleNum, pcrYld=pcrYld)
        totalseqs = 0
        uniqSeqs = 0
        #compute total seq num, unique seq num, and transfer info to x
        for i, seqIdx in enumerate(slctdSeqs):
            uniqSeqs += 1
            totalseqs += int(slctdSeqs[seqIdx][0])
        #initialize matrix to hold info for amplified pool
        x = np.zeros((uniqSeqs, pcrCycleNum+4))
        for i, seqIdx in enumerate(slctdSeqs):
            x[i][0] = seqIdx
            x[i][1] = slctdSeqs[seqIdx][0]
            x[i][2] = slctdSeqs[seqIdx][1]
            x[i][3] = slctdSeqs[seqIdx][2]
        print("number of unique seqs in selected pool prior to amplification: "+str(uniqSeqs))
        print("number of seqs in selected pool prior to amplification: "+str(totalseqs))
        # calculate probabilities of different possible mutation numbers
        mutNumProbs = mut.get_mutation_probabilities_original()
        # compute a discrete distribution of mutation numbers
        mutDist = mut.get_mutation_distribution_original()
        print("Discrete Mutation Distribution has been computed")
    # PCR Amplification
        totalseqs = 0
        # initialize dictionary to keep info on seqs to be mutated
        mutatedPool = {}
        #initialize matrix to hold info for mutation pool
        y = np.zeros((uniqSeqs, seqLength+1))
        # keep track of sequence count after each pcr cycle (except last one)
        seqPop = np.zeros(pcrCycleNum)
        # compute cycle number probabilities for this seq
        cycleNumProbs = np.zeros(pcrCycleNum)
        print("Amplification has started...")
        # for each sequence in the selected pool
        for i, seqIdx in enumerate(list(slctdSeqs)):
            # random PCR with bias using brute force        
            for n in xrange(pcrCycleNum):
                # sequence count after n cycles
                seqPop[n] = x[i][1]
                # amplify count using initial count, polymerase yield, and bias score
                x[i][1] += int(binom(x[i][1], pcrYld+x[i][3]))
            # compute cycle number probabilities
            for s, seqNum in enumerate(seqPop):
                cycleNumProbs[s] = seqNum/np.sum(seqPop)
            # transfer info to x
            for j, cycleNumProb in enumerate(cycleNumProbs):
                x[i][j+4] = cycleNumProb
            # update total num of seqs
            totalseqs += x[i][1]
            # transfer info from x to selection pool
            slctdSeqs[int(x[i][0])] = x[i][1:]
            #tranfer seq index to matrix y
            y[i][0] = x[i][0]
            #if accumulated seq count is greater than 10,000
            if np.sum(seqPop) > 10000:
            # for each possible number of mutations in any seq copy (1-seqLength)
                for mutNum in xrange(seqLength):
                    #approximate the proportion of copies that will be mutated using
                    #corresponding probability p(M=mutNum)
                    y[i][mutNum+1] = mutNumProbs[mutNum+1]*np.sum(seqPop)
            # if seq count is less than 10,000
            else:
                # draw random mutNum from the mutation distribution for each seq copy
                muts = poisson(errorRate*seqLength, int(np.sum(seqPop))) #SLOW STEP
                # remove all drawn numbers equal to zero
                muts = muts[muts != 0]
                # for each non-zero mutation number
                for mutNum in muts:
                    #increment copy number to be mutated
                    y[i][mutNum+1] += 1
        del(x)
        gc.collect()
        print("Amplification carried out")
        print("Sequence selection for mutation has started...")
        #remove all mutation numbers with zero copies to be mutated
        y = y[y[:, 1] != 0]
        #for each seq to be mutated
        for mutInfo in y:
            #add to mutation pool with it's corresponding mutation info
            mutatedPool[int(mutInfo[0])] = mutInfo[1:][mutInfo[1:] != 0]
        del(y)
        gc.collect()
        print("Mutation selection has been carried out")
        print("Mutant generation has started...")
        if(distance == "hamming"):
            # generate mutants and add to the amplfied sequence pool 
            amplfdSeqs = mut.generate_mutants_1D(mutatedPool=mutatedPool, 
                                                 amplfdSeqs=slctdSeqs, 
                                                 aptamerSeqs=aptamerSeqs, 
                                                 alphabetSet=alphabetSet)
            del(slctdSeqs)
            del(mutatedPool)
            gc.collect()
        elif(distance == "basepair"):
            amplfdSeqs = mut.generate_mutants_2D(mutatedPool=mutatedPool,
                                                 amplfdSeqs=slctdSeqs,
                                                 aptamerSeqs=aptamerSeqs,
                                                 alphabetSet=alphabetSet)
            del(slctdSeqs)
            del(mutatedPool)
            gc.collect()
        elif(distance == "loop"):
            amplfdSeqs = mut.generate_mutants_loop(mutatedPool=mutatedPool,
                                                 amplfdSeqs=slctdSeqs,
                                                 aptamerSeqs=aptamerSeqs,
                                                 alphabetSet=alphabetSet)
            del(slctdSeqs)
            del(mutatedPool)
            gc.collect()
        else:
            print("argument given for distance is invalid")
            return
        print("Mutation has been carried out")
        return amplfdSeqs
