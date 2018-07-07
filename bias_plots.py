#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright 2016 Alaa Abdel Latif. All rights reserved.
import sys
import numpy as np
from utils import bias_avg
from scipy.interpolate import spline
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

#This generate a plot of the average bias of total
#and unique sequences during selex
def generate_bias_plot(fileNames, roundNum, seqLength):
    weighted_bias_per_rnd = np.zeros(roundNum)
    bias_per_rnd = np.zeros(roundNum)
    for rnd in xrange(roundNum):
        weighted_bias_per_rnd[rnd], bias_per_rnd[rnd] = bias_avg(fileNames+"_R"+str(rnd+1), 
                                                                 seqLength)
    roundNumAxis = np.linspace(0, roundNum, roundNum)
    plotsList = [weighted_bias_per_rnd, bias_per_rnd]
    fig0, axes = plt.subplots(1, 2)
    basic_colors = ['b', 'g']
    for i, ax in enumerate(axes.reshape(-1)):
        ax.plot(roundNumAxis, plotsList[i], color=basic_colors[i])
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fig0.text(0.5, 0.04, 'Round Number', ha='center')
    fig0.text(0.04, 0.5, 'Average Bias', va='center', rotation='vertical')
    fig0.text(0.3, 0.95, 'Total Sequences', ha='center')
    fig0.text(0.725, 0.95, 'Unique Sequences', ha='center')
    fig0.text(0.07, 0.9, '(a)', ha='center')
    fig0.text(0.507, 0.9, '(b)', ha='center')
    fig0.savefig(str(fileNames)+"_SELEX_Analytics_bias", format='pdf')

#generate_bias_plot('he4_loop_small', 40, 20)

#This generates a plot of the average bias of each affinity group
#during selex
def generate_bias_per_dist_plot(fileNames, roundNum, seqLength, distance):
    weighted_bias_per_rnd = np.zeros((seqLength+1, roundNum))
    bias_per_rnd = np.zeros((seqLength+1, roundNum))
    if(distance=="hamming"):
        for rnd in xrange(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R"+str(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[1:seqLength-1], bias_per_rnd[1:seqLength-1]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0, 0.9, seqLength-1)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in xrange(seqLength-2):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d+1))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist", format='pdf')
        return fig0
    elif(distance=="basepair"):
        for rnd in xrange(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R"+str(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[:seqLength-8], bias_per_rnd[:seqLength-8]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0, 0.9, seqLength-8)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in xrange(seqLength-8):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist", format='pdf')
        return fig0
    elif(distance=="loop"):
        for rnd in xrange(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R"+str(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[:seqLength+1], bias_per_rnd[:seqLength+1]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0.1, 0.9, seqLength+1)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in xrange(seqLength+1):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist", format='pdf')
        return fig0
    else:
        print("Invalid distance metric. Exiting...")
        return 0

#fig = generate_bias_per_dist_plot('he4_loop_small', 40, 20, "loop")
