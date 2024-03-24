#!/bin/bash

#paar, vgrG, etc. This will pull genes that may or may not have extensions
# need to enter a given genome as the first positional parameter.

#pfam05488 PAAR domain
#pfam13665 DUF4150 (PAAR like)

#pfam05954 Phage GPD protein (VgrG analog)
#pfam06715 gp5 repeats (VgrG component)
#pfam04717 bsae V (VgrG component)
#pfam13296 T6SS_VgrG

#TIGR01646 type VI secretion system tip protein VgrG
#TIGR03361 type VI secretion system tip protein TssI/VgrG N-terminal domain

#COG3501 Uncharacterized conserved protein VgrG, implicated in type VI secretion and phage assembly

#pfam05638 HCP
#pfam17642 TssD


grep -H -w 'pfam05954\|pfam05488\|pfam13665\|pfam17642\|pfam05638\|pfam06715\|pfam04717\|pfam13296' /IMG/$1/$1*pfam* >> core_search_outputs/$1_pfam_cores.tsv
grep -H -w 'TIGR01646\|TIGR03361'  /IMG/$1/$1*tigrfam* >> core_search_outputs/$1_tigr_vgrg_cores.tsv
grep -H -w 'COG3501' /IMG/$1/$1*cog* >> core_search_outputs/$1_cog_vgrg_cores.tsv

sed -i 's|.tab.txt:|\t|' core_search_outputs/$1_pfam_cores.tsv
sed -i 's|.tab.txt:|\t|' core_search_outputs/$1_tigr_vgrg_cores.tsv
sed -i 's|.tab.txt:|\t|' core_search_outputs/$1_cog_vgrg_cores.tsv

awk '{print $2}' core_search_outputs/$1_*_cores.tsv | sort -u >> core_search_outputs/$1_core_genes.tsv

seqkit grep -f core_search_outputs/$1_core_genes.tsv /IMG/$1/$1*genes.faa* >> faa_raw/"$1"_t6ss_cores.faa
