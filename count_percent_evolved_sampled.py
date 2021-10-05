# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:25:22 2021

@author: asafl-lab
"""


#%% Sampling by genome
import random
import matplotlib.pyplot as plt


f = open('C:/Users/asafl-lab/Desktop/evolved T6SS/evolved_counts/tree_collapse_groups_0.1')
collapsed = [ eval(line.strip()) for line in f.readlines() ] 


num_samples = 10000

for i in range(num_samples):
    sampled = [ random.sample(group, 1)[0] for group in collapsed ]
    
    sampled_cores = pfam_hits[pfam_hits.genome_id.isin(sampled)].drop_duplicates('gene_oid')
        
    if i == 0:
        evolved_counts = sampled_cores.groupby('core').evolved.value_counts(normalize = True)
    else:
        evolved_counts = pd.concat([evolved_counts, sampled_cores.groupby('core').evolved.value_counts(normalize = True)], axis = 1)
        
evolved_counts_mean = evolved_counts.mean(axis = 1)
evolved_counts_std = evolved_counts.std(axis = 1)




fig, axes = plt.subplots(1,3, figsize= (15, 15))

for i, core in enumerate(cores_list):
        
         axes[i].bar(core, evolved_counts_mean.loc[(core, False)], bottom = evolved_counts_mean.loc[(core, True)], yerr = evolved_counts_std.loc[(core, False)], label = 'Unevolved', alpha=0.85)
         axes[i].bar(core, evolved_counts_mean.loc[(core, True)], yerr = evolved_counts_std.loc[(core, True)], label = 'Evolved', alpha=0.85)   
         if i == 0:
             axes[i].set_ylabel("Percent", size=40)

         axes[i].tick_params(axis='both', labelsize=20)

         
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

# fig.legend(lines, labels)

fig.legend(lines[:2], labels[:2], fontsize = 40, title='Legend', bbox_to_anchor=(1.0001, 1), loc='upper left')

fig.suptitle('Number of evolved VS non-evolved - sampled' + str(num_samples) + ' iterations'  , fontsize=32)


fig.tight_layout()

fig.savefig('sampled by genome_evolved VS nonEvolved counts_'+ str(num_samples)+'_samples', dpi = 150)



#%% Plotting the sampled


cores_list = ['PAAR', 'VGRG', 'HCP']
fig, axes = plt.subplots(1,3, figsize= (15, 15))

for i, core in enumerate(cores_list):
    df_temp = sampled_cores[sampled_cores.core == core].drop_duplicates('gene_oid').evolved.value_counts(normalize=True)
    axes[i].bar(core, df_temp.loc[False], bottom = df_temp.loc[True], label = 'Unevolved', alpha=0.85)
    axes[i].bar(core, df_temp.loc[True], label = 'Evolved', alpha=0.85)
    
    
    
    # axes[i].set_xlabel(core, size=15)
    if i == 0:
        axes[i].set_ylabel("Percent", size=40)
    # axes[i].set_title(core + " in T6SS operons compared to orphan", size = 20)
    # axes[i].legend(loc='upper right', fontsize = 40)
    
    # axes[i].set_xticks(fontsize = 0.001)
    axes[i].tick_params(axis='both', labelsize=20)
    
    
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

# fig.legend(lines, labels)

fig.legend(lines[:2], labels[:2], fontsize = 40, title='Legend', bbox_to_anchor=(1.0001, 1), loc='upper left')


fig.tight_layout()

fig.savefig('evolved VS nonEvolved counts - All', dpi = 150)













