mkdir -p Data_Bases
cd Data_Bases

#downloaded from TAIR
makeblastdb -in TAIR10_cds_20110103_representative_gene_model_updated -out ArabidopsisRepGenes -dbtype nucl

makeblastdb -in /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/B73_LongestTS.fasta -out B73_CDS -dbtype nucl

#removed trailing ".v3.1" from transcript names
makeblastdb -in /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/BTX_623_LongestTS.fasta -out BTX623_CDS -dbtype nucl
makeblastdb -in /panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa -out Mo17_CDS -dbtype nucl

# removed the trailing "m" from transcript gene names in the fasta file
makeblastdb -in /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/phb47.repcds.fa -out PHB47_CDS -dbtype nucl
makeblastdb -in /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa -out PH207_CDS -dbtype nucl
makeblastdb -in /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/W22_LongestTS.fasta -out W22_CDS -dbtype nucl