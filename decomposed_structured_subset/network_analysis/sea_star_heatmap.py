# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 13:14:54 2020

@author: jared
"""

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import json

#fname = "sea_star_Hinf0_medium"
fname = "sea_star_Hinf0_wide_med"

with open(fname + '/sea_star_results.json') as f:
    data = json.load(f)
    
    
cost_primal = np.array(data["cost_primal"])    
cost_primal[np.isnan(cost_primal)] = np.Infinity
time_primal = np.array(data["time_primal"])    
time_dual = np.array(data["time_dual"])    

cost_dual = np.array(data["cost_dual"])
    
RUN_PSD = 0

cones = ["$\mathcal{D}\mathcal{D}$"] + ["$B_{%d}$"%i for i in [1,3,5,8,15, 30, 55]] + ["$\mathbb{S}_+$"]
if RUN_PSD:
    #add psd cone
    #cones = data["cones"]+["psd"]
    cost_primal  = np.append(cost_primal, np.ones([1, cost_primal.shape[1]])*data["cost_psd"], axis = 0)
    time_primal  = np.append(time_primal, np.ones([1, cost_primal.shape[1]])*data["time_psd"], axis = 0)
    
    
    cost_dual  = np.append(cost_dual, np.ones([1, cost_dual.shape[1]])*data["cost_psd"], axis = 0)
    time_dual  = np.append(time_dual, np.ones([1, cost_primal.shape[1]])*data["time_psd"], axis = 0)
    
    
    
    primal_mask = np.zeros_like(cost_primal )
    primal_mask[(cost_primal - data["cost_psd"] >= 1e-4)] = True
    cost_psd = data["cost_psd"]
else:
    #cones = data["cones"]    
    primal_mask = ~np.array(data["sdp_opt_primal"])
    cost_psd = cost_primal[-1, -1]

#plot the results
tick_size = 18
#sup_size = 22
title_size = 22
label_size = 16
cbar_size = 18
annot_size = 16
with sns.axes_style("white"):
    fig = plt.figure(3)
    fig.clf()
    #plt.subplot(1, 2, 1)
    sns.heatmap(time_primal.T/60, mask=primal_mask.T, annot=True, fmt="0.2f", annot_kws={"size":annot_size})
    plt.yticks(np.arange(cost_primal.shape[1])+0.5, data["thresh"], size=tick_size)
    plt.xticks(np.arange(cost_primal.shape[0])+0.5, cones, size=tick_size)
    plt.ylabel('PSD threshold', size = label_size)
    plt.xlabel('Cone Complexity', size = label_size)
    plt.title('Upper bound time (min.) $\gamma = %0.3f$'% cost_psd, size=title_size)
    plt.tight_layout()

    fig2 = plt.figure(4)
    fig2.clf()
    sns.heatmap(cost_dual.T,  cmap="YlGnBu", annot=True, fmt="0.3f", annot_kws={"size":annot_size})
    plt.yticks(np.arange(cost_dual.shape[1])+0.5, data["thresh"], size=tick_size)
    plt.xticks(np.arange(cost_dual.shape[0])+0.5, cones, size=tick_size)
    plt.ylabel('PSD threshold', size = label_size)
    plt.xlabel('Cone Complexity', size = label_size)
    plt.title('Lower bound $\gamma$', size=title_size)

    plt.tight_layout()
    
    fig3 = plt.figure(5)
    fig3.clf()
    sns.heatmap(time_dual.T/60, annot=True, fmt="0.2f", annot_kws={"size":annot_size})
    plt.yticks(np.arange(time_dual.shape[1])+0.5, data["thresh"], size=tick_size)
    plt.xticks(np.arange(time_dual.shape[0])+0.5, cones, size=tick_size)
    plt.ylabel('PSD threshold', size = label_size)
    plt.xlabel('Cone Complexity', size = label_size)
    plt.title('Lower bound time (min.)', size=title_size)

    plt.tight_layout()