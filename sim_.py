#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import sys
import gc
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
from postprocess import dataAnalysis
import utils

#Fetch experiment parameters from the settings file
import ConfigParser
settings = ConfigParser.ConfigParser()
settings.read('settings.init')
#Assign parameters to variables
aptamerType = settings.get('general','selex_type')
aptamerNum = settings.getint('general','aptamer_mode')
aptamerSeq = settings.get('general','reference_aptamer')
seqLength = settings.getint('general','sequence_length')
roundNum = settings.getint('general','number_of_rounds')
outputFileNames = settings.get('general','experiment_name')
samplingSize = settings.getint('general','sampling_size')
post_process = settings.get('general','post_process')
selectionThreshold = settings.getint('selectionparams','scale')
distanceMeasure = settings.get('selectionparams','distance')
stringency = settings.getint('selectionparams','stringency')
pcrCycleNum = settings.getint('amplificationparams','number_of_pcr')
pcrYield = settings.getfloat('amplificationparams','pcr_efficiency')
pcrErrorRate = settings.getfloat('amplificationparams','pcr_error_rate')

# Instantiating classes
Apt = Aptamers()
S = Selection()
Amplify = Amplification()
Mut = Mutation()
# Primary Structure-based Evolution
if(distanceMeasure == "hamming"):
# SELEX simulation based on random aptamer assignment, hamming-based stochastic selection, and
# non-ideal stochastic amplfication with pyrimidine-based bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            del(slctdSeqs)
            gc.collect()
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            del(amplfdSeqs)
            gc.collect()
            print("Selection carried for R"+str(r+1))
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            del(slctdSeqs)
            gc.collect()
            print("Amplification carried for R"+str(r+1))
            print("writing R"+str(r+1)+" seqs to file")
            outFile = outputFileNames + "_R" + str(r+1)
            with open(outFile, 'w') as f:
                for seqIdx in amplfdSeqs:
                    seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                    f.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n')
    print("SELEX completed")

elif(distanceMeasure == "basepair"):
# SELEX simulation based on random aptamer assignment, basepair-based stochastic selection, and
# non-ideal stochastic amplfication with pyrimidine-based bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticBasePairSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")
elif(distanceMeasure == "loop"):
# SELEX simulation based on random aptamer assignment, motif-based stochastic selection, and
# non-ideal stochastic amplfication with pyrimidine-based bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticLoopSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticLoopSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")
else:
    print("Invalid argument for distance measure")
if(post_process):
    print("Data post-processing has started...")
    dataAnalysis(seqLength, roundNum, outputFileNames, post_process, distanceMeasure)
    print("Data post-processing is complete")
    print("The Simulation has ended")
else:
    print("The Simulation has ended without post-processing")
