#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 16:26:27 2023

@author: raharinirina
"""

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import pdb
import pickle
import seaborn as sns
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import os

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u

#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)
    
def Display_Stacked(t, t_dates, Y, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 32, save_to = "test", xval = "x", yval = "f(x)", 
            linewidth = 3, palette = None, linestyle = None, color = None, ax = None, fig = None, linecolor = None ,alpha = 0.5, mode = None, get_file = False):
    
    if fig == None:
        PreFig(xsize = xysize[0], ysize = xysize[1])
        fig = plt.figure(figsize = figsize)
        if ax == None:
            ax = fig.add_subplot(1, 1, 1)
    
    if linestyle is None:
        linestyle = ["-"]*Y.shape[0]
    
    
    Y = np.concatenate((np.zeros((1,Y.shape[1])), Y), axis = 0)
    
    Y = np.cumsum(Y, axis = 0)
    
    if mode == "freq":
        Y = 100*np.divide(Y, (Y[-1, :][np.newaxis, :]), out = np.zeros(Y.shape), where = Y[-1, :][np.newaxis, :]!=0)

    elif mode == "prop":
        Y = np.divide(Y, (Y[-1, :][np.newaxis, :]), out = np.zeros(Y.shape), where = Y[-1, :][np.newaxis, :]!=0)
    
    
    labels = [" "] + list(labels) 
    if palette is None:
        for i in range(1, Y.shape[0]):
            if color is None:
                ax.fill_between(t, Y[i-1, :], Y[i, :], label = labels[i], alpha = alpha)
            else:
                ax.fill_between(t, Y[i-1, :], Y[i, :], label = labels[i], color = color[i-1], alpha = alpha)
            
            if linecolor is not None:
                plt.plot(t, Y[i, :], color = linecolor, linewidth = linewidth)

    else:
        col = sns.color_palette(palette, Y.shape[0])
        for i in range(1, Y.shape[0]):
            ax.fill_between(t, Y[i-1, :], Y[i, :], label = labels[i], color = col[i-1], alpha = alpha)
        
            if linecolor is not None:
                plt.plot(t, Y[i, :], color = linecolor, linewidth = linewidth)
    
    if is_log:
        plt.ylabel("%s (log$_{10}$)"%yval, fontsize = labsize)
    else:
        plt.ylabel("%s"%yval, fontsize = labsize)  
    
    plt.xlabel(xval, fontsize = labsize)
    plt.legend(loc = (1.2, 0), fontsize = labsize, ncols = np.ceil(len(labels)/20).astype(int))
    
    if len(t)>200:
        pp = 7*4
    else:
        pp = 14
        
    t_ticks = np.array(t)[::pp]
    t_ticks_labels = np.array(t_dates)[::pp]
    if t[-1] not in t_ticks:
        t_ticks = np.append(t_ticks, t[-1])
        t_ticks_labels = np.append(t_ticks_labels, t_dates[-1])
    
    ax.set_xticks(t_ticks)
    ax.set_xticklabels(t_ticks_labels, rotation = 45, horizontalalignment = "right")
        
    ax.set_xlim((t[0], t[-1]))


    if save_to[-3:] == "pdf":
        ### save figure in pdf ###
        pdf = PdfPages(save_to)
        pdf.savefig(fig, bbox_inches = "tight")
        if not get_file:
            pdf.close()
        fig.savefig(save_to[:-4]+".svg", bbox_inches = "tight")
        plt.close()

    else:
        plt.savefig(save_to)
        pdf = None
        
    if get_file:
        return fig, ax, pdf
    else:
        return fig, ax


def Display_Envelops(t, t_dates, Y, Z, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 32, save_to = "test", xval = "x", yval = "f(x)", 
            linewidth = 3, palette = None, linestyle = None, color = None, ax = None, fig = None, linecolor = None ,alpha = 0.5, mode = None, get_file = False, yfmt = None):
    
    if fig == None:
        PreFig(xsize = xysize[0], ysize = xysize[1])
        fig = plt.figure(figsize = figsize)
        if ax == None:
            ax = fig.add_subplot(1, 1, 1)
    
    if linestyle is None:
        linestyle = ["-"]*Y.shape[0]
    
    labels = list(labels) 
    if palette is None:
        for i in range(Y.shape[0]):
            if color is None:
                ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], alpha = alpha)
            else:
                if i == 0:
                    ax.plot(t, Y[i, :], label = labels[i], color = color[i], linewidth = 5, alpha = 0.4)
                    ax.plot(t, Z[i, :], color = color[i], linewidth = 5, alpha = 0.4)
                else:
                    ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], color = color[i], alpha = alpha)

    else:
        col = list(sns.color_palette(palette, 2*(Y.shape[0]-1))[::2]) + ["#1f77b4"]
        for i in range(Y.shape[0]):
            if i == 0:
                ax.plot(t, Y[i, :], label = labels[i], color = col[i], linewidth = 5, alpha = 0.4)
                ax.plot(t, Z[i, :], color = col[i], linewidth = 5, alpha = 0.4)
            else:
                ax.fill_between(t, Y[i, :], Z[i, :], label = labels[i], color = col[i], alpha = alpha)
        
    
    if is_log:
        plt.ylabel("%s (log$_{10}$)"%yval, fontsize = labsize)
    else:
        plt.ylabel("%s"%yval, fontsize = labsize)  
    
    plt.xlabel(xval, fontsize = labsize)
    #plt.legend(loc = (1.2, 0), fontsize = labsize, ncols = np.ceil(len(labels)/15).astype(int))
    plt.legend(loc = (0.25, -1), fontsize = labsize, ncols = 2)
    if len(t)>200:
        pp = 7*4
    else:
        pp = 14
        
    t_ticks = np.array(t)[::pp]
    t_ticks_labels = np.array(t_dates)[::pp]
    if t[-1] not in t_ticks:
        t_ticks = np.append(t_ticks, t[-1])
        t_ticks_labels = np.append(t_ticks_labels, t_dates[-1])
    
    ax.set_xticks(t_ticks)
    ax.set_xticklabels(t_ticks_labels, rotation = 45, horizontalalignment = "right")
    ax.set_xlim((t[0], t[-1]))
    ax.set_ylim((40e6, 120e6))
    if yfmt is not None:
        y_ticks = np.arange(40e6, 120e6, 10*10**yfmt)
        if np.abs(np.max(Y) - y_ticks[-1]) > 10**yfmt:
            y_ticks = np.append(y_ticks, (np.max(Y)//10**yfmt)*10**yfmt)
            
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(["%d"%(y_ticks[i]//10**yfmt) for i in range(0, len(y_ticks)-1)]+["%d x $ 10^{6}$"%(y_ticks[-1]//10**yfmt)])
        #ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p:format(int(x)//10, ",")))
    
    if save_to[-3:] == "pdf":
        ### save figure in pdf ###
        pdf = PdfPages(save_to)
        pdf.savefig(fig, bbox_inches = "tight")
        if not get_file:
            pdf.close()
        fig.savefig(save_to[:-4]+".svg", bbox_inches = "tight")
        plt.close()

    else:
        plt.savefig(save_to)
        pdf = None
        
    if get_file:
        return fig, ax, pdf
    else:
        return fig, ax
    
    
def Display(t, t_dates, Y, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 32, save_to = "test", xval = "x", yval = "f(x)", 
            linewidth = 3, palette = None, linestyle = None, color = None, ax = None, fig = None, linecolor = None ,alpha = 0.5, mode = None, get_file = False, yfmt = None):
    
    if fig == None:
        PreFig(xsize = xysize[0], ysize = xysize[1])
        fig = plt.figure(figsize = figsize)
        if ax == None:
            ax = fig.add_subplot(1, 1, 1)
    
    if linestyle is None:
        linestyle = ["-"]*Y.shape[0]
    
    labels = list(labels) 
    if palette is None:
        for i in range(Y.shape[0]):
            if color is None:
                ax.plot(t, Y[i, :], label = labels[i], alpha = alpha, linewidth = linewidth, linestyle = linestyle[i])
            else:
                ax.plot(t, Y[i, :], label = labels[i], linewidth = linewidth, color = color[i], alpha = alpha, linestyle = linestyle[i])

    else:
        col = list(sns.color_palette(palette, 2*(Y.shape[0]))[::2])
        for i in range(Y.shape[0]):
            ax.plot(t, Y[i, :], label = labels[i],linewidth = linewidth, color = col[i], alpha = alpha, linestyle = linestyle[i])
        
    
    if is_log:
        plt.ylabel("%s (log$_{10}$)"%yval, fontsize = labsize)
    else:
        plt.ylabel("%s"%yval, fontsize = labsize)  
    
    plt.xlabel(xval, fontsize = labsize)
    plt.legend(loc = (1.2, 0), fontsize = labsize, ncols = np.ceil(len(labels)/15).astype(int))
    
    if len(t)>200:
        pp = 7*4
    else:
        pp = 14
        
    t_ticks = np.array(t)[::pp]
    t_ticks_labels = np.array(t_dates)[::pp]
    if t[-1] not in t_ticks:
        t_ticks = np.append(t_ticks, t[-1])
        t_ticks_labels = np.append(t_ticks_labels, t_dates[-1])
    
    ax.set_xticks(t_ticks)
    ax.set_xticklabels(t_ticks_labels, rotation = 45, horizontalalignment = "right")
    ax.set_xlim((t[0], t[-1]))
    ax.set_ylim((-0.3e6, 6e6))
    
    if yfmt is not None:
        #y_ticks = np.arange(0, np.max(Y), 10*10**yfmt)
        y_ticks = np.arange(0, np.ceil(np.max(Y)), 10*10**yfmt)
        #if np.abs(np.max(Y) - y_ticks[-1]) > 10**yfmt:
            #y_ticks = np.append(y_ticks, (np.max(Y)//10**yfmt)*10**yfmt)
        y_ticks = np.append(y_ticks, (np.ceil(np.max(Y)/10**yfmt))*10**yfmt)
            
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([0]+["%d"%(y_ticks[i]//10**yfmt) for i in range(1, len(y_ticks)-1)]+["%d x $ 10^{6}$"%(y_ticks[-1]//10**yfmt)])
        #ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p:format(int(x)//10, ",")))
    
    if save_to[-3:] == "pdf":
        ### save figure in pdf ###
        pdf = PdfPages(save_to)
        pdf.savefig(fig, bbox_inches = "tight")
        if not get_file:
            pdf.close()
        fig.savefig(save_to[:-4]+".svg", bbox_inches = "tight")
        plt.close()

    else:
        plt.savefig(save_to)
        pdf = None
        
    if get_file:
        return fig, ax, pdf
    else:
        return fig, ax


col1 = "coolwarm"
col2 = "coolwarm"
"Mean susceptible from infection Landscape: Without vaccination"
file1, file2 = str(sys.argv[1]), str(sys.argv[2])
S_mean_df = pd.read_csv(file1)
S_all_mean_orig = S_mean_df.to_numpy()[:, (S_mean_df.columns != "Days")&(S_mean_df.columns != "Unnamed: 0")].astype(float)
t_dates_orig = S_mean_df["Days"].tolist()

"Mean susceptible from infection + vaccination Landscape: With vaccination"
S_mean_df_2 = pd.read_csv(file2)
S_all_mean_orig_2 = S_mean_df_2.to_numpy()[:, (S_mean_df_2.columns != "Days")&(S_mean_df_2.columns != "Unnamed: 0")].astype(float)
t_dates_orig_2 = S_mean_df_2["Days"].tolist()


"Vaccination data "
vacc_infos = pd.read_csv(sys.argv[3])
vacc_names = vacc_infos.columns[(vacc_infos.columns != "date")&(vacc_infos.columns != "Unnamed: 0")].tolist()
Counts = vacc_infos.to_numpy()[:, (vacc_infos.columns != "date")&(vacc_infos.columns != "Unnamed: 0")].astype(float)
weights = np.divide(Counts, np.sum(Counts, axis = 1)[:, np.newaxis], out = np.zeros(Counts.shape), where = np.sum(Counts, axis = 1)[:, np.newaxis]!= 0)
t_vacc = vacc_infos["date"].tolist() ### the same for all the vacc_name and equale to date column in vacc_infos


#S_vacc_dir = sys.argv[4] ## IL Vacc vs Spikes
#S_all_dir = sys.argv[5] ## IL Spike vs Vacc
file1 = open(sys.argv[6], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()

N_pop = float(sys.argv[7])

date_start = sys.argv[8]
date_end = sys.argv[9]

if date_start not in t_dates_orig:
    date_start = t_dates_orig[0]
if date_end not in t_dates_orig:
    date_end = t_dates_orig[-1]
t_dates = np.array(t_dates_orig)#[t_dates_orig.index(date_start):t_dates_orig.index(date_end)+1]
S_all_mean = S_all_mean_orig#[t_dates_orig.index(date_start):t_dates_orig.index(date_end)+1, :]

d_min = list(t_dates).index(date_start)
d_max = list(t_dates).index(date_end) + 1


# plotting proportions
PreFig(xsize = 40, ysize = 40)
fig = plt.figure(figsize = (15, 7))
ax = fig.add_subplot(1, 1, 1) 
filename = sys.argv[-1]+"/Vacc_counts_ver1.pdf"


Vacc_Counts_ver1 = Counts.T

fig, ax = Display(t_vacc, t_vacc, np.cumsum(Vacc_Counts_ver1, axis = 1), is_log = False, labels = vacc_names, color = ["black", "grey"], figsize = None, save_to = filename, xval = "dates", yval = "Counts (CumSum)", 
                        linewidth = 8, linestyle = ["-", "--"], ax = ax, fig = fig, alpha = 1, mode = None, yfmt = 6)

cdic = {"dates":t_vacc}
cdic.update({vacc_names[i]:np.cumsum(Vacc_Counts_ver1, axis=1)[i, :] for i in range(len(vacc_names))})


c_df = pd.DataFrame(cdic)
c_df.to_csv(sys.argv[-1]+"/Vacc_counts_ver1.csv")

### plot run_ver
S2 = np.zeros(S_all_mean.shape)
S_mask = np.zeros(S_all_mean.shape).astype(bool)
for i in range(len(t_dates)):
    if t_dates[i] in t_dates_orig_2:
        w_i = list(t_dates_orig_2).index(t_dates[i])
        S2[i, :] = S_all_mean_orig_2[w_i, :]
    else:
        S_mask[i, :] = True

S2 = ma.masked_array(S2, mask = S_mask)
S_mask_sum = np.all(S_mask, axis = 1)
S2_min, S2_max = np.min(S2, axis = 1), np.max(S2, axis = 1)

Y_2_min = []
Y_2_max = []

Y_2_min.append(S2_min)
Y_2_max.append(S2_max)

Y_2_min.append(np.min(S_all_mean, axis = 1))
Y_2_max.append(np.max(S_all_mean, axis = 1))

Y_2_min = np.array(Y_2_min)
Y_2_max = np.array(Y_2_max)

Y_2_min = ma.masked_array(Y_2_min, mask = np.row_stack(tuple([S_mask_sum for i in range(Y_2_min.shape[0])])))
Y_2_max = ma.masked_array(Y_2_max, mask = np.row_stack(tuple([S_mask_sum for i in range(Y_2_min.shape[0])])))

# plotting
PreFig(xsize = 40, ysize = 40)
fig = plt.figure(figsize = (15, 7))
ax = fig.add_subplot(1, 1, 1) 

filename = sys.argv[-1]+"/Susceptible_Trends_ALL_vs_Vacc.pdf"

t = np.arange(len(t_dates))
if str(sys.argv[-1][-4:]) == "ver1":
    run_ver = "ALL vs VACC only"
elif str(sys.argv[-1][-4:]) == "ver2":
    run_ver = "With vaccination"
fig, ax = Display_Envelops(t[d_min:d_max], t_dates[d_min:d_max], Y_2_min[:, d_min:d_max], Y_2_max[:, d_min:d_max], is_log = False, labels = [run_ver,  "Without vaccination"], color = ["red", "green", "#1f77b4"], figsize = None, save_to = filename, xval = "dates", yval = "Mean $\mathbb{E}$[Susceptible]", 
                          linewidth = 3, ax = ax, fig = fig,alpha = 0.3, mode = None, yfmt = 6)

