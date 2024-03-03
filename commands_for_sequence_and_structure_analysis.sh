# Commands used for sequence as well as structural clustering, and for structure-based search

#---------------

#Create an mmseqs database of all fasta amino acid files
mmseqs createdb *.faa all_cores_db

# Full length clustering with mmseqs
mmseqs cluster  --cov-mode 0 -c 0.9  --min-seq-id 0.5 all_cores_db all_cores_cluster_0.9COV_0.5_ID_db tmp

#Make into a readable tsv
mmseqs createtsv all_cores_db all_cores_db all_cores_cluster_0.9COV_0.5_ID_db all_cores_cluster_0.9COV_0.5_ID_db.tsv

#---------------

#Create an mmseqs db for only the C-termini
mmseqs createdb c_termini_minusSpiralsminusPF.faa c_termini_minusSpiralsminusPF_DB

#Cluster C-termini by sequence
mmseqs cluster  --cov-mode 0 -c 0.85  --min-seq-id 0.5 c_termini_minusSpiralsminusPF_DB c_termini_minusSpiralsminusPF_DB_Clustered_covMode0_cov0.85_id0.5 tmp

#Make into a readable tsv
mmseqs createtsv  c_termini_minusSpiralsminusPF_DB  c_termini_minusSpiralsminusPF_DB c_termini_minusSpiralsminusPF_DB_Clustered_covMode0_cov0.85_id0.5 c_termini_minusSpiralsminusPF_DB_Clustered_covMode0_cov0.85_id0.5.tsv

#Filter for C-termini representatives faa
cat c_termini_minusSpiralsminusPF_DB_Clustered_covMode0_cov0.85_id0.5.tsv | awk '{print $1}' | sort -u > c_termini_minusSpiralsminusPF_representative_IDs.tsv
seqkit grep -f c_termini_minusSpiralsminusPF_representative_IDs.tsv c_termini_minusSpiralsminusPF.faa > c_termini_minusSpiralsminusPF_representatives.faa

#---------------
#pfam search

hmmscan --domtblout c_termini_alt.domtblout Pfam-A.hmm ./c_termini_minusSpiralsminusPF_representatives.faa

#---------------

#Predicting 3D structure of proteins (/ folding proteins with immunities)
colabfold_batch --num-recycle 20 --recycle-early-stop-tolerance=0.5 --stop-at-score 90 --model-type alphafold2_multimer_v3 --rank multimer c_termini_minusSpiralsminusPF_representatives.faa output"

#---------------

#create foldseek structure database from pdb files predicted by alphafold
foldseek createdb ./c_termini_minusSpiralsminusPF_representatives_rank1s_structure_predicitons/*.pdb c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB

#cluster by structure using foldseek
foldseek cluster --cov-mode 0 -c 0.85 -e 1e-6 \
c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB \
c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB_FOLDSEEK_CLUSTER_covmode0_cov0.85_e1e-6_DB \
tmp

#Make into a readable tsv
foldseek createtsv c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB \
c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB \
c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB_FOLDSEEK_CLUSTER_covmode0_cov0.85_e1e-6_DB \
c_termini_minusSpiralsminusPF_representatives_rank1s_FOLDSEEK_DB_FOLDSEEK_CLUSTER_covmode0_cov0.85_e1e-6_DB.tsv


