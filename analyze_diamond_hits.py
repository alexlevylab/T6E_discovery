import numpy as np
import pandas as pd

hits_and_surrounding_gff = pd.read_csv("round1_diamond_hits_plus_upstream_and_downstream_labeledGroups_cleaned.gff", sep = '\t', 
usecols = [0, 1, 4, 5, 7, 9, 10], header = None, names = ['genome_ID', 'scaffold', 'start_coord', 'end_coord', 'strand', 'gene_oid', 'hit_group'])

hits_and_surrounding_gff['rank_in_group'] = hits_and_surrounding_gff.groupby('hit_group')['start_coord'].rank()

###

diamond_hits_raw = pd.read_csv("round1_total__vs__t6ss_db4_DIAMOND_subject_geneoids_unique", sep = '\t', header = None, names = ['gene_oid'])

diamond_hits_raw['raw_diamond_hit'] = True

diamond_hits_raw['gene_oid'] = diamond_hits_raw['gene_oid'].astype(int)

#print(diamond_hits_raw)
###

mer1 = hits_and_surrounding_gff.merge(diamond_hits_raw, how = 'left', on = 'gene_oid')

mer1['raw_diamond_hit'].fillna(False, inplace = True)

#print(mer1)
###
scaffolds_of_interest = mer1.loc[mer1['raw_diamond_hit'] == True][['gene_oid', 'scaffold', 'hit_group']].drop_duplicates()

scaffolds_of_interest.rename(columns = {"scaffold" : "scaffold_of_interest"}, inplace = True)

#print(scaffolds_of_interest)

mer1 = mer1.merge(scaffolds_of_interest, on = ['hit_group', 'gene_oid'], how = 'left')

mer1['scaffold_of_interest'] = mer1.groupby('hit_group')['scaffold_of_interest'].fillna(method = 'ffill')
mer1['scaffold_of_interest'] = mer1.groupby('hit_group')['scaffold_of_interest'].fillna(method = 'bfill')

#remove edge cases by removing any that are at the end of a scaffold
mer1 = mer1.groupby('hit_group').filter(lambda x : (x['scaffold'] == x['scaffold_of_interest']).all())

print(mer1.shape)

###

original_core_gene_hits = pd.read_csv("t6ss_clusters_with_gt_50_percent_of_genomes_with_all3_TssJLM_ALL_VgrG_PAAR_HCP.pfam.half_cleaned2", sep = '\t',
header = None)

original_core_gene_hits = original_core_gene_hits[[1, 9, 10]]

original_core_gene_hits.rename(columns = {1: 'gene_oid', 9: 'core_pfam', 10: 'core_name'}, inplace = True)

original_core_gene_hits['t6ss_core_in_gene'] = True

original_core_gene_hits['gene_oid'] = original_core_gene_hits['gene_oid'].astype(int)

mer2 = mer1.merge(original_core_gene_hits, on = 'gene_oid', how = 'left')

mer2['t6ss_core_in_gene'].fillna(False, inplace = True)


#################
# classification of diamond hits:
# 1) total orphans -> diamond hits with no core in it, nor in surroundings
# 2) next to cores -> diamond hits with no core, but yes has a core next to it
# 3) has a core -> diamond hit with a core ... exclude

#1
total_orphans = mer2.groupby('hit_group').filter(lambda x: (x['t6ss_core_in_gene'] == False).all())

print("diamond hits that are true orphans", total_orphans['hit_group'].nunique())


#2
diamond_hit_next_to_core = mer2.groupby('hit_group').filter(lambda x: ((x['t6ss_core_in_gene'] == False) & (x['raw_diamond_hit'] == True)).any() & ((x['t6ss_core_in_gene'] == True) & (x['raw_diamond_hit'] == False)).any() & ~((x['t6ss_core_in_gene'] == True) & (x['raw_diamond_hit'] == True)).any())

#3
diamond_hit_has_core = mer2.groupby('hit_group').filter(lambda x: ((x['t6ss_core_in_gene'] == True) & (x['raw_diamond_hit'] == True)).any())

#######

total_orphans = total_orphans.loc[total_orphans['raw_diamond_hit'] == True]
diamond_hit_next_to_core = diamond_hit_next_to_core.loc[diamond_hit_next_to_core['raw_diamond_hit'] == True]
diamond_hit_has_core = diamond_hit_has_core.loc[diamond_hit_has_core['raw_diamond_hit'] == True]

total_orphans['total_orphans'] = True

diamond_hit_next_to_core['diamond_hit_next_to_core'] = True

diamond_hit_has_core['diamond_hit_has_core'] = True

one_two_combined = pd.concat([total_orphans, diamond_hit_next_to_core, diamond_hit_has_core])

#######

# import cluster IDs and query IDs

members = pd.read_csv("round2_65_ID_80_cover_members.csv")

members.rename(columns = {"gene_oid" : "query", "count" : "c_term_cluster_group"}, inplace = True)

diamond_query_and_subject = pd.read_csv("round1_total__vs__t6ss_db4_DIAMOND", sep = '\t', usecols = [0, 1], header = None, names = ['query', 'subject'])

diamond_reference = diamond_query_and_subject.merge(members, how = 'left', on = 'query')

######

final = diamond_reference.merge(one_two_combined, how = 'left', left_on = 'subject', right_on = 'gene_oid')

final2 = final.loc[~final['gene_oid'].isnull()]

print(final2.groupby("c_term_cluster_group")['diamond_hit_has_core'].value_counts())
print(final2.groupby("c_term_cluster_group")['diamond_hit_next_to_core'].value_counts())
print(final2.groupby("c_term_cluster_group")['total_orphans'].value_counts())



final2.to_csv("diamond_results_processed.csv", index = False)