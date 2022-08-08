# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 13:14:54 2020

@author: jared
"""

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
#import pandas as pd

#from matplotlib import rc
#import matplotlib.pylab as plt

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

#uncons_cost = pd.read_csv('../LR120_uncons_dual_cost.csv', sep = ',', header=None)
#uncons_cost = np.genfromtxt('../LR120_uncons_dual_cost.csv', delimiter = ',')
#uncons_time = np.genfromtxt('../LR120_uncons_dual_time.csv', delimiter = ',')

uncons_cost = np.genfromtxt('./LR_output_uncons_v2_cost.csv', delimiter = ',')
uncons_time = np.genfromtxt('./LR_output_uncons_v2_time.csv', delimiter = ',')


cons_cost = np.genfromtxt('../LR120_output_box_1_2_cost.csv', delimiter = ',')
cons_time = np.genfromtxt('../LR120_output_box_1_2_time.csv', delimiter = ',')

uncons_cost[np.isnan(uncons_cost)] = np.inf

uncons_mask = np.zeros_like(uncons_time)
#cost_gap = abs(uncons_cost-uncons_cost[-1, -1])
#cost_gap_mask = np.isnan(cost_gap) | (cost_gap >= 1e-3)
uncons_mask[abs(uncons_cost-uncons_cost[-1, -1]) >= 1e-3] = True
uncons_mask[np.isnan(uncons_cost)] = True

cons_mask = np.zeros_like(cons_time)
#cons_mask[abs(cons_cost-cons_cost[-1, -1]) >= 1e-3] = True
cons_mask[0:4, 0:3] = True

#uniform_data = np.random.rand(15, 20)
#
#mask = np.zeros_like(uniform_data)
#mask[np.triu_indices_from(mask)] = True
#

dd = '$\mathcal{D}\mathcal{D}$'
sdd = '$\mathcal{S}\mathcal{D}\mathcal{D}$'
psd = '$\mathbb{S}_+$'

#uncons_inner = ['$B_{%d}$' % i for i in ( 2, 3, 5, 6, 11, 20, 30, 40)]
#
#uncons_cone = [dd, sdd] +  uncons_inner + [psd]
#
##uncons_cone = (dd, sdd, 2, 3, 5, 6, 11, 20, 30, 40, psd)
#uncons_thresh = (0, 12, 45, 100)
#
#cons_inner = ['$B_{%d}$' % i for i in (2, 3, 5,	6,	10,	15,	20)]
#cons_cone = [dd, sdd] +  cons_inner + [psd]
#
##cons_cone = (dd, sdd, 2, 3, 5,	6,	10,	15,	20,	psd)
#cons_thresh = (0, 5, 12, 45, 100)

cone_inner = ['$B_{%d}$' % i for i in (2, 3, 5,	6,10,15,	20, 30, 40)]

cone = [dd, sdd] + cone_inner + [psd]
thresh = (0, 5, 12, 45, 100)



#uncons_note = np.zeros_like(uncons_mask)
#ind_uncons_max =  np.unravel_index(np.argmax(uncons_time - 1e9*uncons_mask, axis=None), uncons_time.shape)
#ind_uncons_min =  np.unravel_index(np.argmin(uncons_time + 1e9*uncons_mask, axis=None), uncons_time.shape)
#
#
#uncons_note[ind_uncons_max] = True
#uncons_note[ind_uncons_min] = True
#
#
#
#cons_note = np.zeros_like(cons_mask)
#
#ind_cons_max =  np.unravel_index(np.argmax(cons_time - 1e9*cons_mask, axis=None), cons_time.shape)
#ind_cons_min =  np.unravel_index(np.argmin(cons_time + 1e9*cons_mask, axis=None), cons_time.shape)
##np.argmin(cons_time + 1e9*cons_mask)
#
#cons_note[ind_cons_max] = True
#cons_note[ind_cons_min] = True

tick_size = 16
#sup_size = 22
title_size = 22
label_size = 16
cbar_size = 18
annot_size = 16

#%%
with sns.axes_style("white"):
#    fig = plt.figure(3)
#    fig.clf()
#    plt.subplot(1, 2, 1)
#    sns.heatmap(uncons_time.T, mask=uncons_mask.T)
#    plt.yticks(np.arange(uncons_cost.shape[1])+0.5, thresh, size=tick_size)
#    plt.xticks(np.arange(uncons_cost.shape[0])+0.5, cone, size=tick_size)
#    plt.ylabel('PSD threshold', size = label_size)
#    plt.xlabel('Cone Complexity', size = label_size)
#    plt.title('Unconstrained $f^*=-110.2$', size=title_size)
#
#    
#    plt.subplot(1, 2, 2)
#    sns.heatmap(cons_time.T, mask=cons_mask.T )
#    plt.yticks(np.arange(cons_cost.shape[1])+0.5, thresh, size=tick_size)
#    plt.xticks(np.arange(cons_cost.shape[0])+0.5, cone, size=tick_size)
#    plt.ylabel('PSD threshold', size = label_size)
#    plt.xlabel('Cone Complexity', size = label_size)
#    plt.title('Constrained $f^*=4939.1$', size=title_size)
#    
#    
##    fig.suptitle('Time (seconds) to find Lower Bounds' , size=sup_size)
#    #sns.heatmap(uniform_data, mask=mask)
#    plt.tight_layout(w_pad = 25)
    
    fig = plt.figure(3)
    fig.clf()
    #plt.subplot(1, 2, 1)
    sns.heatmap(uncons_time.T/60, mask=uncons_mask.T, annot=True, fmt="0.2f",annot_kws={"size":annot_size})
    plt.yticks(np.arange(uncons_cost.shape[1])+0.5, thresh, size=tick_size)
    plt.xticks(np.arange(uncons_cost.shape[0])+0.5, cone, size=tick_size)
    plt.ylabel('PSD threshold', size = label_size)
    plt.xlabel('Cone Complexity', size = label_size)
    plt.title('Unconstrained $f^*=-110.2$ time (min.)', size=title_size)

    
    #plt.subplot(1, 2, 2)
    fig = plt.figure(4)
    fig.clf()
    sns.heatmap(cons_time.T/60, mask=cons_mask.T, annot=True, fmt="0.2f",annot_kws={"size":annot_size})
    plt.yticks(np.arange(cons_cost.shape[1])+0.5, thresh, size=tick_size)
    plt.xticks(np.arange(cons_cost.shape[0])+0.5, cone, size=tick_size)
    plt.ylabel('PSD threshold', size = label_size)
    plt.xlabel('Cone Complexity', size = label_size)
    plt.title('Constrained $f^*=4939.1$ time (min.)', size=title_size)
    
    
#    fig.suptitle('Time (seconds) to find Lower Bounds' , size=sup_size)
    #sns.heatmap(uniform_data, mask=mask)
    plt.tight_layout(w_pad = 25)

plt.show()