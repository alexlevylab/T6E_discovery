# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:37:50 2021

@author: asafl-lab
"""
#%% PLOTTING ORPHANS

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize = (80, 30))
for idx, core_name in enumerate(core_dict):
    
    df = pfam_hits[pfam_hits.core.isin([core_name])].drop_duplicates('gene_oid', keep = 'first')
    
    
    df_evolved_around = df[df.evolved == True].num_T6SS_around.to_list()
    df_nonEvolved_around = df[df.evolved != True].num_T6SS_around.to_list()
    
    # bins = [0,1] + list(range(2,23,2))
    bins = range(0,23, 1)
    
    ax[idx].hist(df_nonEvolved_around, bins=bins, alpha=0.5, label="Non Evolved " + core_name, density=True)
    ax[idx].hist(df_evolved_around, bins=bins, alpha=0.5, label="Evolved " + core_name, density=True)
    
    
    ax[idx].set_xlabel("Number of T6SS components in " + core_name +" neighbourhood", size=50)
    ax[idx].set_ylabel("Count", size=50)
    ax[idx].set_title(core_name + " in T6SS operons compared to orphan", size = 65)
    ax[idx].legend(loc='upper right', fontsize = 40)
    
    # ax[idx].set_xticks(bins, fontsize = 30)
    ax[idx].tick_params(axis='both', labelsize=30)
    
    
# fig.savefig('all cores histogram evolved histogram_fixed VGRG_density.png', dpi = 150)