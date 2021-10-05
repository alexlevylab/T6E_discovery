# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 12:12:12 2021

@author: asafl-lab
"""

import random
import pandas as pd
import matplotlib.pyplot as plt


plt.style.use('seaborn')

pfam_hits = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/evolved domains pie charts/pfam_hits_data_with_cores_for_graphs.tsv')

f = open('C:/Users/asafl-lab/Desktop/evolved T6SS/evolved_counts/tree_collapse_groups_0.1')
collapsed = [ eval(line.strip()) for line in f.readlines() ] 

core_dict = {'VGRG' : ['Gp5_C', 'Phage_GPD', 'Phage_base_V', 'T6SS_Vgr', 'DUF2345'], 'PAAR' : ['DUF4150', 'PAAR_motif'], 'HCP' : ['T6SS_HCP'] }

num_samples = 10




fig, axes = plt.subplots(1,3, figsize= (100, 100))
c = 0

for i, core in enumerate(core_dict):
    for n in range(num_samples):
        sampled = [ random.sample(group, 1)[0] for group in collapsed ]
        
        sampled_cores = pfam_hits[pfam_hits.genome_id.isin(sampled)]
        

        df = sampled_cores[sampled_cores.core == core]
        
        df_evolved_core = df[(df.evolved == True)]
        evolved_domains = set(df[~df.pfam_name.isin(core_dict[core])].pfam_name.unique())
        evolved_with_unknown_domain = (df_evolved_core.groupby('gene_oid').pfam_name.unique().apply(lambda x: set(x)).apply(lambda x: len(x.intersection(evolved_domains))) == 0).value_counts().loc[True]
        
        # 
        df = df[~df.pfam_name.isin(core_dict[core])]
        
        
        if n == 0:
            counts = df.pfam_name.value_counts()
            counts['Unknown'] = evolved_with_unknown_domain
        
        else:
            
            counts2 = df.pfam_name.value_counts()
            counts2['Unknown'] = evolved_with_unknown_domain
            
            counts = pd.concat([counts, counts2], axis = 1)
    
    
    counts = counts.mean(axis = 1).sort_values(ascending = False)
    
    explode = [ 0 if c != 'Unknown' else 0.1 for c in counts.index ]
    
    axes.flat[i].pie(counts, textprops={'fontsize': 40} , explode = explode)#wedgeprops = {'edgecolor': 'white'})
    axes.flat[i].tick_params(axis='x', labelrotation= 90, labelsize = 25)
    axes.flat[i].set_title(core, fontdict = {'fontsize': 200})
    if i!=1:
        axes.flat[i].legend(counts.index,title="evolved domains",loc="upper right", prop={'size': 50}, bbox_to_anchor= (1.3,1))
    else:
        axes.flat[i].legend(counts.index,title="evolved domains",loc="upper right", prop={'size': 50}, bbox_to_anchor= (1.3,1.6))
    
    
    fig.suptitle('Evolved pfam domains per core - sampled' + str(num_samples) + ' iterations'  , fontsize=32)


    fig.tight_layout()
    
    # fig.savefig('Evolved pfam domains per core - sampled '+ str(num_samples), dpi = 150)


    
    

    
    
    