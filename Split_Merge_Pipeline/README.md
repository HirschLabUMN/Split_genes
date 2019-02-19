# Split/merge gene identification pipeline
There are are two seperate pipelines that were used in this analysis  
* Split/merge identification without anchors using nucmer.  
* Split/merge identification using an anchor based pipleline (Dagchainer).
## Preprocessing of the data
* Here is the list of metadata files that were used for this analysis and thier locations:  
```
samples:
  b73: {bed: last/b73-phb47/b73-phb47.q.bed, cds: /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.cds.all.fa,
    fasta: /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa,
    gff: /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.gff3,
    prot: /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.pep.all.fa,
    rep_bed: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/b73-ph207/b73.repcds.bed,
    rep_cds: /home/hirschc1/broha006/software/synmap/b73-phb47/b73.repcds.fa, rep_fasta: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/b73-ph207/b73.repcds.fa}
  mo17: {cds: /home/maize/shared/databases/genomes/Zea_mays/GCA_003185045.1_Zm-Mo17-REFERENCE-CAU-1.0_cds_from_genomic_AB.fna,
    fasta: /home/maize/shared/databases/genomes/Zea_mays/mo17_10chr_gbrowse_1August2018.fa,
    gff: /home/maize/shared/databases/genomes/Zea_mays/GCA_003185045.1_Zm-Mo17-REFERENCE-CAU-1.0_genomic.groomed.gff,
    prot: /home/maize/shared/databases/genomes/Zea_mays/GCA_003185045.1_Zm-Mo17-REFERENCE-CAU-1.0_protein_AB.faa,
    rep_bed: /panfs/roc/groups/14/hirschc1/broha006/software/synmap/w22-mo17/mo17.repcds.bed,
    rep_fasta: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa}
  ph207: {cds: /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa,
    fasta: /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa,
    gff: /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.gff3,
    prot: /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.protein.fa,
    rep_bed: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/ph207-w22/ph207.repcds.bed,
    rep_fasta: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/ph207-phb47/ph207.repcds.fa}
  phb47: {cds: /home/maize/shared/databases/genomes/Zea_mays/PHB47/ZmaysvarPHB47v1.1.allTrs.cds.fa,
    fasta: /home/maize/shared/databases/genomes/Zea_mays/PHB47/Zea_mays_var_PHB47.mainGenome.fasta,
    gff: /home/maize/shared/databases/genomes/Zea_mays/PHB47/Zmaysvar.PHB47v1.1.gene_exons.gff3,
    mRNA: /home/maize/shared/databases/genomes/Zea_mays/PHB47/ZmaysvarPHB47v1.1.allTrs.fa,
    prot: /home/maize/shared/databases/genomes/Zea_mays/PHB47/ZmaysvarPHB47v1.1.allTrs.pep.fa,
    rep_bed: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-phb47/phb47.repcds.bed,
    rep_fasta: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/b73-phb47/phb47.repcds.fa}
  w22: {bed: foo.txt, cds: /home/maize/shared/databases/genomes/Zea_mays/W22/zea_maysw22.coding.fasta,
    fasta: /home/maize/shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta,
    gff: /home/maize/shared/databases/genomes/Zea_mays/W22/zea_maysw22_core_fixchr.gff,
    prot: /home/maize/shared/databases/genomes/Zea_mays/W22/zea_maysw22.protein.fasta,
    rep_bed: /panfs/roc/groups/14/hirschc1/broha006/software/synmap/w22-mo17/w22.repcds.bed,
    rep_fasta: /panfs/roc/groups/13/stuparr/mich0391/software/synmap/w22-b73/w22.repcds.fa}
```
One of the first things that we needed to do was process the files to pull out the longest representative transcript.  
* PH207 only has once transcript per gene so we did not need to modify that file

This is the list of the longest transcript gff and exon files that I created/used for makeing a blastable database and downstream analysis  
```
#Longest transcript gff files
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Mo17_Longest_transcript_Exons.gff
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Zea_mays.AGPv4.33_Longest_transcript_Exons.gff
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/zea_maysw22_core_fixchr_Longest_Transcript_Exon.gff
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/PHB47_Longest_TS.gff
/home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.gff3

# representative cds files
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/B73_LongestTS.fasta
/home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/W22_LongestTS.fasta
/panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/BTX_623_LongestTS.fasta
/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/phb47.repcds.fa
```
* If you need to extract the longest transcript from a file, use the Parse_GFF_LongestTranscript.py script the make a new fasta file with the representative CDS using Get_L_fasta.py
```python3
python3 Preprocess/Parse_GFF_LongestTranscript.py -g /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.gff3 \
    -o Zea_mays.AGPv4.33_Longest_transcript_Exons.gff \
    -i Parent

# from the Zea_mays.AGPv4.33_Longest_transcript_Exons.gff output file, you will need to subset transcript names into its own file format
# You can find a quick python example script on how to to this here
    # Keep in mind that you will need to modify the split field on the 9th column to subset the right transcript name
python3 Preprocess/get_longest_TSname.py > B73_longest_TS.txt

python3 Preprocess/Get_L_fasta.py \
    -i B73_longest_TS.txt \
    -f Zea_mays.B73_RefGen_v4.cds.all.fa \
    -o B73_LongestTS.fasta
```

* Once the longest representive transcript and CDS fasta has been extracted, you can make a blastable database
```bash
# here are the cammands I used to make all of the blastable databases
bash Preprocess/MakeBlastDB.sh
```

## Nucmer based pipeline to identify split/merged
* The first step is to run nucmer on all chromosomal fasta files
```bash
# all by all propogator for number commands
bash Nucmer_All_by_All/mummer_propogator.sh > all_mummer_commands.txt
# using a series pipes and greps, I broke the commands up into smaller files for running on MSI
# For example
bash Nucmer_All_by_All/mummer_propogator.sh | grep "p B73_" > B73_nucmer_commands.sh
# Script to run it on MSI from mesabi node
qsub Nucmer_All_by_All/Mummer_Run_B73.sh

# once all of the mummer commands have finished, we use a script to reformat the headers on output file
# The command will run on all files within the directory the command is run
bash Nucmer_All_by_All/mummer_post_process.sh
```
* Once all databases and nucmer output files have run, we then run the all by all pipeline
```python3
# modify the propogator file with your respective files and paths
# The "SPLIT" variable is the character used to split gene number with chromosome so we can see if genes are numerically next to each other

# Example command
python All_by_All_Blast.py \
    -q /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/B73_genes.gff3 \
    -s /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/PH207_genes_noV1.gff3 \
    -b /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/PH207_CDS \
    -l /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/B73_LongestTS.fasta \
    -o B73_PH207_AllbyAll_res.txt \
    -n /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Mummer_Results/B73_PH207_c1000.fil.coords \
    -g a

# Command propogator for all line combinations
bash Nucmer_All_by_All/all_by_all_propogator.sh > Nucmer_All_by_All/aba_commands.txt

# The python All_by_All_Blast.py command can be run simutanious in parallel
# If you wanted to run three commands on MSI you can run a command like
module load parallel
parallel --jobs 3 < Nucmer_All_by_All/aba_commands.txt
```