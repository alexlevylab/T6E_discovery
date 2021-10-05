# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 15:34:57 2021

@author: asafl-lab
"""
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
# import scipy.stats as stats




col_names = ['gene_oid',
 'length_aa',
 'query_start',
 'query_end',
 'subj_start',
 'subj_end',
 'evalue',
 'bit_score',
 'pfam_id',
 'pfam_name',
 'domain_length']




pfam_hits = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/data/all_pfam_hits_exec_No_RVTDomain_wDUF4150.tsv', names = col_names)
pfam_hits['genome_id'] = pfam_hits.gene_oid.str.split('/', expand=True)[8]
pfam_hits.gene_oid = pfam_hits.gene_oid.str.replace('.pfam.tab.txt:', '/').str.split('/', expand=True)[10]
pfam_hits.gene_oid = pfam_hits.gene_oid.astype('int64')

#Removing ends of scaffolds
genes_in_end_of_scaffolds = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/end_of_scaffolds/end_of_scaffolds.txt')
pfam_hits = pfam_hits[~pfam_hits.gene_oid.isin(genes_in_end_of_scaffolds.gene_oid)]

#Removing phage genes
phage_hits = pd.read_csv('C:/Users/asafl-lab/Desktop/evolved T6SS/removing_phages/all_phage_hits.tsv', sep=' ')
pfam_hits = pfam_hits[~pfam_hits.gene_oid.isin(phage_hits[phage_hits.hits != 0].gene_oid)]



core = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/data/all_HCP_PAAR_VGRG_hit_wPF13665.tsv', names = col_names)
core['genome_id'] = core.gene_oid.str.split('/', expand=True)[8]
core.gene_oid = core.gene_oid.str.replace('.pfam.tab.txt:', '/').str.split('/', expand=True)[10]
core.gene_oid = core.gene_oid.astype('int64')

# adding DUF2345 as core domain
core = pd.concat([core, pfam_hits[pfam_hits.pfam_name == 'DUF2345']])


# add COG data for VGRG
core_COG_data = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/vgrg/add_COGs/vgrg_core_COG_data.tsv')
core = pd.concat([core, core_COG_data]).sort_values(['gene_oid', 'query_end'])
core.length_aa = core.groupby('gene_oid').ffill().length_aa


# add TIGRFAM for vgrg
core_TIGRFAM_data = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/vgrg/add_TIGRFAM/vgrg_core_TIGRFAM_data.tsv')
core = pd.concat([core, core_TIGRFAM_data]).sort_values(['gene_oid', 'query_end'])
core.length_aa = core.groupby('gene_oid').ffill().length_aa


core['evolved'] = core.gene_oid.map(core.groupby('gene_oid').last()['length_aa'] - core.groupby('gene_oid').last()['query_end'] > 70)



last_core_domain_coord = core.sort_values(['gene_oid', 'query_end']).groupby('gene_oid').last().query_end




#remove COG rows
core = core[~core.pfam_id.str.contains('COG')]
core = core[~core.pfam_id.str.contains('TIGR')]




core['core_length'] = core.subj_end - core.subj_start
core_length = core[['pfam_name', 'domain_length']].drop_duplicates().drop_duplicates('pfam_name')
core_length['length_threshold'] = core_length.domain_length * 0.6
core = pd.merge(core, core_length[['pfam_name', 'length_threshold']], how = 'left', on = 'pfam_name')

partial_domain = core[core.core_length < core.length_threshold].gene_oid





pfam_hits = pfam_hits[~pfam_hits.gene_oid.isin(partial_domain)]
pfam_hits = pd.merge(pfam_hits, core[['gene_oid', 'evolved']], how = 'left', on = 'gene_oid')

core_dict = {'VGRG' : ['Gp5_C', 'Phage_GPD', 'Phage_base_V', 'T6SS_Vgr', 'DUF2345'], 'PAAR' : ['DUF4150', 'PAAR_motif'], 'HCP' : ['T6SS_HCP'] }
core_domains =  ['Gp5_C', 'Phage_GPD', 'Phage_base_V', 'T6SS_Vgr', 'DUF2345', 'DUF4150', 'PAAR_motif', 'T6SS_HCP']



pfam_hits['core'] = ''
for core_type in core_dict:
    genes_w_core = pfam_hits[pfam_hits.pfam_name.isin(core_dict[core_type])].gene_oid
    pfam_hits.loc[pfam_hits.gene_oid.isin(genes_w_core), 'core'] = core_type



#Removing small VGRGs
pfam_hits = pfam_hits[~pfam_hits.gene_oid.isin(pfam_hits.loc[(pfam_hits.core == 'VGRG') & (pfam_hits.length_aa < 500)].gene_oid)]


neighbourhood = pd.read_table('C:/Users/asafl-lab/Desktop/evolved T6SS/neighbourhood_orphan_or_operon.csv')
neighbourhood.loc[neighbourhood.num_T6SS_around >= 5, 'neighbourhood'] = 'operon'

pfam_hits = pd.merge(pfam_hits, neighbourhood, how = 'left', on = 'gene_oid').sort_values(['gene_oid', 'query_end'])
pfam_hits['domain_archi'] = pfam_hits.gene_oid.map(pfam_hits.groupby('gene_oid').pfam_name.unique().str.join('~'))


# removing genes that dont start with core domain
pfam_hits = pfam_hits[~pfam_hits.gene_oid.isin(pfam_hits.groupby('gene_oid').first()[~pfam_hits.groupby('gene_oid').first().pfam_name.isin(core_domains)].index)]



# setting evolved VGRG with GMM - THRESHHOLD FOR EVOLVD: mean of unevolved gaussian + 1 std (which is covariance**0.5) ###
vgrg = pfam_hits[pfam_hits.core == 'VGRG']
X = np.array(vgrg.drop_duplicates('gene_oid').length_aa.to_list()).reshape(-1, len(vgrg.drop_duplicates('gene_oid').length_aa.to_list())).T

gm = GaussianMixture(n_components=2, random_state=0, covariance_type = 'full').fit(X)

means_gaussians =  gm.means_
covariances_gaussians = gm.covariances_
weights__gaussians =  gm.weights_

threshhold_length_for_evolved = round((means_gaussians[0][0] + covariances_gaussians[0][0]**0.5)[0]) #setting threshold for evolved vgrg




pfam_hits['last_core_domain_coord'] = pfam_hits.gene_oid.map(last_core_domain_coord)


pfam_hits.loc[(pfam_hits.core == 'VGRG') & (pfam_hits.length_aa > threshhold_length_for_evolved) & (pfam_hits.length_aa - pfam_hits.last_core_domain_coord > 70), 'evolved'] = True
pfam_hits.loc[(pfam_hits.core == 'VGRG') & (pfam_hits.length_aa <= threshhold_length_for_evolved), 'evolved'] = False



# pfam_hits.to_csv('pfam_hits_data_with_cores_for_graphs.tsv', index = False, sep = '\t')


# pfam_hits.sort_values(['gene_oid', 'query_start']).to_csv('C:/Users/asafl-lab/Desktop/evolved T6SS/data/all_pfam_hits_exec_No_RVTDomain_wDUF4150_removed_scaffold_ends_noPhage_removedSmallVGRG_smaller_than_500_Evolved_VGRG_GMM.tsv', index=False, sep='\t')





#%% Adding fasta and slicing evolved domains


from Bio import SeqIO

# all_cores = []

# for seq in SeqIO.parse('C:/Users/asafl-lab/Desktop/evolved T6SS/all_cores.fa', 'fasta'):
#     if int(seq.id) in pfam_hits.gene_oid.unique():
#         all_cores.append(seq)




evolved = pfam_hits[pfam_hits.evolved == True].sort_values(['gene_oid', 'query_end'])


evolved_coords = evolved[['gene_oid', 'core', 'last_core_domain_coord']].drop_duplicates().set_index('gene_oid')


vgrg_sliced_seqs = []
hcp_sliced_seqs = []
paar_sliced_seqs = []


all_sliced_evolved_seqs = []

evolved_coords_l = list(evolved_coords.index)

c = 0 
for seqrecord in SeqIO.parse('C:/Users/asafl-lab/Desktop/evolved T6SS/evolved_cores.fa', 'fasta'):
    if int(seqrecord.id) in evolved_coords_l:
        core, domain_end = evolved_coords.loc[int(seqrecord.id)].core, evolved_coords.loc[int(seqrecord.id)].last_core_domain_coord
        if (core == 'VGRG') and (domain_end < 797):
            domain_end = 797
        seqrecord.seq = seqrecord.seq[domain_end:]
        
        if core == 'VGRG':
            vgrg_sliced_seqs.append(seqrecord)
        elif core == 'PAAR':
            paar_sliced_seqs.append(seqrecord)
        else:
            hcp_sliced_seqs.append(seqrecord)

        all_sliced_evolved_seqs.append(seqrecord)
    
        if c%1000 == 0:
            print(c)
        c+=1

    
SeqIO.write(vgrg_sliced_seqs, 'evolved_vgrg.fasta', 'fasta')
SeqIO.write(hcp_sliced_seqs, 'hcp_vgrg.fasta', 'fasta')
SeqIO.write(paar_sliced_seqs, 'paar_vgrg.fasta', 'fasta')
SeqIO.write(all_sliced_evolved_seqs, 'all_vgrg.fasta', 'fasta')






























