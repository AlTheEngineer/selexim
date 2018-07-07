#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import sys
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import cm
import Distance

D = Distance.Distance()

#This generate the main plots from the simulation results
#The plots include changes in total and unique sequence numbers, in average distance
#and changes in the average distance of each affinity group
def dataAnalysis(seqLength, roundNum, outputFileNames, plots, distanceMeasure, aptSeq=None, aptStruct=None, aptLoop=None):
    avgDist_per_rnd = np.zeros(roundNum)
    weighted_avgDist_per_rnd = np.zeros(roundNum)
    total_seqs_freqs = np.zeros(roundNum)
    uniq_seqs_freqs = np.zeros(roundNum)
    distFreqs = np.zeros((roundNum, seqLength+5))
    weighted_distFreqs = np.zeros((roundNum, seqLength+5))
    for rnd in xrange(roundNum):
        total_seq_num = 0
        uniq_seq_num = 0
        distance = 0
        weighted_distance = 0
        with open(outputFileNames+"_R"+str(rnd+1)) as SELEX_round:
            for line in SELEX_round:
                columns = line.split()
                distance += int(columns[2])
                weighted_distance += int(columns[1])*int(columns[2])
                total_seq_num += int(columns[1])
                uniq_seq_num += 1
                distFreqs[rnd][int(columns[2])] += 1
                weighted_distFreqs[rnd][int(columns[2])] += int(columns[1])
        avgDist_per_rnd[rnd] = int(distance/uniq_seq_num)
        weighted_avgDist_per_rnd[rnd] = int(weighted_distance/total_seq_num)
        total_seqs_freqs[rnd] = total_seq_num
        uniq_seqs_freqs[rnd] = uniq_seq_num
    for rnd in xrange(roundNum):
	    for i in xrange(seqLength+5):
		    distFreqs[rnd][i] /= uniq_seqs_freqs[rnd]
		    weighted_distFreqs[rnd][i] /= total_seqs_freqs[rnd]
    with open(outputFileNames+"_processed_results", 'w') as p:
        for rnd in xrange(roundNum): 
            p.write(str(int(total_seqs_freqs[rnd]))+'\t')
            p.write(str(int(uniq_seqs_freqs[rnd]))+'\t')
            p.write(str(int(avgDist_per_rnd[rnd]))+'\t')
            p.write(str(int(weighted_avgDist_per_rnd[rnd]))+'\t')
            for l in xrange(seqLength+1):
                p.write(str(int(distFreqs[rnd][l]))+'\t')
                p.write(str(int(weighted_distFreqs[rnd][l]))+'\t')
            p.write('\n')
   # If the user requested generating plots
    if(plots==True):
        # If Hamming distances were used
        if(distanceMeasure=="hamming"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            fig0, axes = plt.subplots(2, 2)
            colormap = plt.cm.gist_ncar
            plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
            colors = [colormap(i) for i in np.linspace(0, 0.9, seqLength-1)]
            basic_colors = ['b', 'g', 'r', 'y']
            for i, ax in enumerate(axes.reshape(-1)):
                ax.plot(roundNumAxis, plotsList[i], color=basic_colors[i])
                if(i==0 or i==1):
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            fig0.text(0.5, 0.04, 'Round Number', ha='center')
            fig0.text(0.04, 0.25, 'Average Distance', va='center', rotation='vertical')
            fig0.text(0.04, 0.725, 'Frequency', va='center', rotation='vertical')
            fig0.text(0.3, 0.95, 'Total Sequences', ha='center')
            fig0.text(0.725, 0.95, 'Unique Sequences', ha='center')
            fig0.text(0.07, 0.9, '(a)', ha='center')
            fig0.text(0.07, 0.475, '(b)', ha='center')
            fig0.text(0.507, 0.9, '(c)', ha='center')
            fig0.text(0.507, 0.475, '(d)', ha='center')
            fig0.savefig(str(outputFileNames)+"_SELEX_Analytics_distance", format='pdf')
            fig1, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                if(i==1):
                    for d in range(2):
                        ax.plot(roundNumAxis, distFreqs[:,d+(3*i)+1], 
                                label='d = '+str(d+(3*i)+1))
                elif(i==2):
                    for d in range(4):
                        ax.plot(roundNumAxis, distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                else:
                    for d in range(3):
                        ax.plot(roundNumAxis, distFreqs[:,d+(3*i)+1], 
                                label='d = '+str(d+(3*i)+1))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*3+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig1.text(0.5, 0.95, 'Unique Sequences', ha='center')
            fig1.text(0.09, 0.9, '(a)', ha='center')
            fig1.text(0.365, 0.9, '(b)', ha='center')
            fig1.text(0.64, 0.9, '(c)', ha='center')
            fig1.text(0.09, 0.475, '(d)', ha='center')
            fig1.text(0.365, 0.475, '(e)', ha='center')
            fig1.text(0.64, 0.475, '(f)', ha='center')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs", format='pdf')
            # weighted fractional sequency plots
            fig2, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                if(i==1):
                    for d in range(2):
                        ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)+1], 
                                label='d = '+str(d+(3*i)+1))
                elif(i==2):
                    for d in range(4):
                        ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                else:
                    for d in range(3):
                        ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)+1], 
                                label='d = '+str(d+(3*i)+1))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*3+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig2.text(0.5, 0.04, 'Round Number', ha='center')
            fig2.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig2.text(0.5, 0.95, 'Total Sequences', ha='center')
            fig2.text(0.09, 0.9, '(a)', ha='center')
            fig2.text(0.365, 0.9, '(b)', ha='center')
            fig2.text(0.64, 0.9, '(c)', ha='center')
            fig2.text(0.09, 0.475, '(d)', ha='center')
            fig2.text(0.365, 0.475, '(e)', ha='center')
            fig2.text(0.64, 0.475, '(f)', ha='center')
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs", format='pdf')
        # If Base Pair distances were used
        elif(distanceMeasure=="basepair"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            # Plots for Average Distances and Sequence Frequencies 
            colormap = plt.cm.gist_ncar
            plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
            colors = [colormap(i) for i in np.linspace(0, 0.9, seqLength-7)]
            fig0, axes = plt.subplots(2, 2)
            basic_colors = ['b', 'g', 'r', 'y']
            for i, ax in enumerate(axes.reshape(-1)):
                ax.plot(roundNumAxis, plotsList[i], color=basic_colors[i])
                if(i==0 or i==1):
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            fig0.text(0.5, 0.04, 'Round Number', ha='center')
            fig0.text(0.04, 0.25, 'Average Distance', va='center', rotation='vertical')
            fig0.text(0.04, 0.725, 'Frequency', va='center', rotation='vertical')
            fig0.text(0.3, 0.95, 'Total Sequences', ha='center')
            fig0.text(0.725, 0.95, 'Unique Sequences', ha='center')
            fig0.text(0.07, 0.9, '(a)', ha='center')
            fig0.text(0.07, 0.475, '(b)', ha='center')
            fig0.text(0.507, 0.9, '(c)', ha='center')
            fig0.text(0.507, 0.475, '(d)', ha='center')
            fig0.savefig(str(outputFileNames)+"_SELEX_Analytics_distance", format='pdf')
            # Plots for distance analytics
            fig1, axes = plt.subplots(3, 2)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(2):
                    ax.plot(roundNumAxis, distFreqs[:,d+(2*i)], label='d = '+str(d+(2*i)))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*2+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.get_yaxis().get_offset_text().set_x(0.05)
                ax.get_yaxis().get_offset_text().set_size(7)
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig1.text(0.5, 0.95, 'Unique Sequences', ha='center')
            fig1.text(0.09, 0.9, '(a)', ha='center')
            fig1.text(0.507, 0.9, '(b)', ha='center')
            fig1.text(0.09, 0.62, '(c)', ha='center')
            fig1.text(0.507, 0.62, '(d)', ha='center')
            fig1.text(0.09, 0.35, '(e)', ha='center')
            fig1.text(0.507, 0.35, '(f)', ha='center')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs", format='pdf')
            # Plot for distance analytics weighted by frequencies
            fig2, axes = plt.subplots(3, 2)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(2):
                    ax.plot(roundNumAxis, weighted_distFreqs[:,d+(2*i)], label='d = '+str(d+(2*i)))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*2+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.get_yaxis().get_offset_text().set_x(0.05)
                ax.get_yaxis().get_offset_text().set_size(7)
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig2.text(0.5, 0.04, 'Round Number', ha='center')
            fig2.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig2.text(0.5, 0.95, 'Total Sequences', ha='center')
            fig2.text(0.09, 0.9, '(a)', ha='center')
            fig2.text(0.507, 0.9, '(b)', ha='center')
            fig2.text(0.09, 0.62, '(c)', ha='center')
            fig2.text(0.507, 0.62, '(d)', ha='center')
            fig2.text(0.09, 0.35, '(e)', ha='center')
            fig2.text(0.507, 0.35, '(f)', ha='center')
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs", format='pdf')
        elif(distanceMeasure=="loop"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            fig0, axes = plt.subplots(2, 2)
            colormap = plt.cm.gist_ncar
            plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
            colors = [colormap(i) for i in np.linspace(0, 0.9, seqLength+3)]
            basic_colors = ['b', 'g', 'r', 'y']
            for i, ax in enumerate(axes.reshape(-1)):
                ax.plot(roundNumAxis, plotsList[i], color=basic_colors[i])
                if(i==0 or i==1):
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            fig0.text(0.5, 0.04, 'Round Number', ha='center')
            fig0.text(0.04, 0.25, 'Average Distance', va='center', rotation='vertical')
            fig0.text(0.04, 0.725, 'Frequency', va='center', rotation='vertical')
            fig0.text(0.3, 0.95, 'Total Sequences', ha='center')
            fig0.text(0.725, 0.95, 'Unique Sequences', ha='center')
            fig0.text(0.07, 0.9, '(a)', ha='center')
            fig0.text(0.07, 0.475, '(b)', ha='center')
            fig0.text(0.507, 0.9, '(c)', ha='center')
            fig0.text(0.507, 0.475, '(d)', ha='center')
            fig0.savefig(str(outputFileNames)+"_SELEX_Analytics_distance", format='pdf')
            fig1, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                if(i==5):
                    for d in range(6):
                        ax.plot(roundNumAxis, distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                else:
                    for d in range(3):
                        ax.plot(roundNumAxis, distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*3+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig1.text(0.5, 0.95, 'Unique Sequences', ha='center')
            fig1.text(0.09, 0.9, '(a)', ha='center')
            fig1.text(0.365, 0.9, '(b)', ha='center')
            fig1.text(0.64, 0.9, '(c)', ha='center')
            fig1.text(0.09, 0.475, '(d)', ha='center')
            fig1.text(0.365, 0.475, '(e)', ha='center')
            fig1.text(0.64, 0.475, '(f)', ha='center')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs", format='pdf')
            # weighted fractional sequency plots
            fig2, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                if(i==5):
                    for d in range(6):
                        ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                else:
                    for d in range(3):
                        ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)], 
                                label='d = '+str(d+(3*i)))
                for j, line in enumerate(ax.lines):
                    line.set_color(colors[i*3+j])
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.tick_params(axis='x', labelsize=5)
                ax.tick_params(axis='y', labelsize=5)
                ax.legend(prop={'size':6})
            fig2.text(0.5, 0.04, 'Round Number', ha='center')
            fig2.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig2.text(0.5, 0.95, 'Total Sequences', ha='center')
            fig2.text(0.09, 0.9, '(a)', ha='center')
            fig2.text(0.365, 0.9, '(b)', ha='center')
            fig2.text(0.64, 0.9, '(c)', ha='center')
            fig2.text(0.09, 0.475, '(d)', ha='center')
            fig2.text(0.365, 0.475, '(e)', ha='center')
            fig2.text(0.64, 0.475, '(f)', ha='center')
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs", format='pdf')
            if(aptSeq != None and aptStruct != None and aptLoop != None):
                dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
                uniq_dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
                for rnd in xrange(roundNum):
                    with open(outputFileNames+"_R"+str(rnd+1)+"_dist_components_results", 'r') as p:
                        for line in p:
                            row = line.split()
                            loop_dist = float(row[0])
                            bp_dist = float(row[1])
                            count = float(row[2])
                            uniq = float(row[3])
                            dist_matrx[rnd][int(loop_dist)][int(bp_dist)] += int(count)
                            uniq_dist_matrx[rnd][int(loop_dist)][int(bp_dist)] += int(uniq)
                fig3, axes = plt.subplots(2, 4)
                for i, ax in enumerate(axes.reshape(-1)):
                    if(i == 0):
                        cax = ax.imshow(dist_matrx[i+1], interpolation='nearest', 
                                                                cmap=cm.coolwarm)
                        cbar = fig3.colorbar(cax, ticks=[np.min(uniq_dist_matrx[i+1]),np.max(uniq_dist_matrx[i+1])], ax=ax)
                        ax.set_title('Round '+str(i+1))
                    else:
                        cax = ax.imshow(dist_matrx[i*5], interpolation='nearest', 
                                                                cmap=cm.coolwarm)
                        cbar = fig3.colorbar(cax, ticks=[np.min(uniq_dist_matrx[i*5]),np.max(uniq_dist_matrx[i*5])], ax=ax)
                        ax.set_title('Round '+str(i*5))
                fig3.savefig(str(outputFileNames)+"_SELEX_Analytics_dist_heatmap", format='pdf')
                fig3.text(0.5, 0.98, 'Total Sequences', ha='center')

        else:
            return
#TEST
#fig = dataAnalysis(20, 40, "he4_loop_small", True, "loop", "GTACGACAGTCATCCTACAC", "(((.((......)).)))..", "CAGTCA")
#PUT IN SEPARATE POST-PROCESSING SCRIPT
'''
if(distanceMeasure == 'loop'):
        if(aptSeq != None and aptStruct != None and aptLoop != None):
            dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
            uniq_dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
            for rnd in xrange(roundNum):
                with open(outputFileNames+"_R"+str(rnd+1)) as SELEX_round:
                    for line in SELEX_round:
                        row = line.split()
                        seq = str(row[0])
                        count = int(row[1])
                        loop_dist, bp_dist = D.loop_components_func(aptSeq, aptStruct, 
                                                                  aptLoop, seq, 
                                                                  seqLength)
                        dist_matrx[rnd][loop_dist][bp_dist] += count
                        uniq_dist_matrx[rnd][loop_dist][bp_dist] += 1
                with open(outputFileNames+"_R"+str(rnd+1)+"_dist_components_results", 'w') as p:
                        for loop_dist in xrange(seqLength+1):
                            for bp_dist in xrange(seqLength-6):
                                p.write(str(loop_dist)+'\t'+str(bp_dist)+'\t'+str(dist_matrx[rnd][loop_dist][bp_dist])+'\t'+str(uniq_dist_matrx[rnd][loop_dist][bp_dist])+'\n')
    else:
        print("The sequence, structure, and loop region of the reference aptamer were not given. Thus, heat maps of the Loop-based and BP distances will not be generated.")
'''
