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
* Once all databases and nucmer output files have run, we then run the all by all pipeline (14 ish hours to run on one core)
* This pipeline is set to use a 500kb flanking window (lines 357-360 in script)
* This pipeline will also prioritize by e-value then size of blast match if there is a tie
* You must also be able to call blasnt from the command line
* -g is the gene delimiter to split the ZMxyz from the numerical order to see if genes are within "5" of each other
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

## Comparing reciprocal split merge files
* I did not make a script for this, but I made a series of Python functions that you could either wrap in a script or run in iPython
* Within Split_merge_compare.py you would run these functions interactily on the output files from the All_By_All_Blast.py files in both directions
* using 4 cores takes about 15 seconds
```python
#Load up all of these functions
import pandas as pd
import multiprocessing
from functools import partial
import numpy as np



def one_to_one_compare(comparison_1, comparison_2):
    """Compare one to one gene key files from an
       all by all blast comparsion"""
    total = comparison_1.shape[0]
    interm = 0
    one2one = 0
    for row in comparison_1.itertuples():
        if row.gene_type == "one_to_one_mapping":
            interm += 1
            Q_gene = row.Query_gene
            S_gene = row.Sytentic_genes.split(",")[0]
            # get the index of the query gene in second file using subject gene var
            idx = comparison_2[comparison_2.Query_gene.isin([S_gene])].index.tolist()
            # check to see if the index is empty
            if idx:
                if comparison_2.at[idx[0], "gene_type"] == "one_to_one_mapping":
                    comp_2_S_gene = comparison_2.at[idx[0], "Sytentic_genes"].split(",")[0]
                    if comp_2_S_gene == Q_gene:
                        one2one += 1
    return(total, interm, one2one)


# Split gene iteration function
def recip_compare(s_file, c_file):
    """Compare split and merge gene calls across reciprocal blast files"""
    recip_list = []
    cross_match = 0
    for row in s_file.itertuples():
        Q_gene = row.Query_gene
        adj_gene_interm = row.adjacent_genes.split(";")
        # list comp to strip metadata and keep gene name
        adj_genes = [i.split(",")[0] for i in adj_gene_interm]
        counter = 0
        for adj_gene in adj_genes:
            Q_search = c_file[c_file.Query_gene.isin([adj_gene])].index.tolist()
            # check to to see that the is one to one match first
            if Q_search and c_file.at[Q_search[0], "gene_type"] =="one_to_one_mapping":
                # Q_search will always results in a list size of 1
                S_search = c_file.at[Q_search[0], "Sytentic_genes"].split(",")[0]
                if S_search == Q_gene:
                    counter += 1
        if counter == len(adj_genes):
            cross_match += 1
            recip_list.append(row)
    recip_df = pd.DataFrame(recip_list)
    return(recip_df)

# Create a function to open up the file and set up parallelization
def All_By_All_compare(i_file1, i_file2, cores):
    """Parallel calls for one to one all by all blast, remember we are using
    cores and not threads"""
    num_processes = int(cores)
    compare_1 = pd.read_csv(i_file1, sep="\t", index_col=False)
    compare_2 = pd.read_csv(i_file2, sep="\t", index_col=False)
    chunks = np.array_split(compare_1, num_processes)
    # pool.map will only take one arg so set up partial fill
    parallel = partial(one_to_one_compare, comparison_2=compare_2)
    pool = multiprocessing.Pool(processes=num_processes)
    result_list = pool.map(parallel, chunks)
    # our fuction returns lists and we want each item to sum accord to pos
    result = [sum(i) for i in zip(*result_list)]
    return(result)


# check confirmation of split vs merged genes in both direction
def split_compare(i_file1, i_file2, cores):
    """Parallel calls for split merge compare from by all blast in reciprocal
    direcitons, remember we are using cores and not threads"""
    num_processes = int(cores)
    file1 = pd.read_csv(i_file1, sep="\t", index_col=False)
    file2 = pd.read_csv(i_file2, sep="\t", index_col=False)
    f1_split = file1[file1.gene_type.isin(["adjacent_genes_syntenic"])]
    f2_split = file2[file2.gene_type.isin(["adjacent_genes_syntenic"])]

    # set up for parallelization
    chunks_1 = np.array_split(f1_split, num_processes)
    # pool.map will only take one arg so set up partial fill
    parallel_1 = partial(recip_compare, c_file=file2)
    pool = multiprocessing.Pool(processes=num_processes)
    # Will return list of pd dfs, need to concat
    recip_1 = pd.concat(pool.map(parallel_1, chunks_1))

    # set up for parallelization
    chunks_2 = np.array_split(f2_split, num_processes)
    parallel_2 = partial(recip_compare, c_file=file1)
    recip_2 = pd.concat(pool.map(parallel_2, chunks_2))

    # merge from both tested directions
    merged_recip = pd.concat([recip_1, recip_2])
    return(merged_recip)


# exmpale run commands
# This command will return the total number of genes that are one to one
# to each other in both directions
JMM_ABA_B73_Mo17 = All_By_All_compare("B73_Mo17_AllbyAll_res.txt",
                                      "Mo17_B73_AllbyAll_res.txt",
                                      4)

# This command will return a dataframe with reciprocal split/merged in agreement
JMM_SC_B73_Mo17 = split_compare("B73_Mo17_AllbyAll_res.txt",
                                "Mo17_B73_AllbyAll_res.txt",
                                4)

# You would then write this to a file using a command such as this (leave index in for next script)
JMM_SC_B73_Mo17.to_csv("500kb_B73_Mo17_recip_SM.txt", sep="\t")
```
* Then you can manually scan the files for any interesting candidates.

## Arabidopsis CDS comparison
* To see which genes are represented as split or merged in Arabidopsis run the "QueryDB_CDS.py" script
* I propgated a file called arabidopsis_blast_cmd.txt containing all of the previous commands
* FULL DISCLOSURE: The -o command does nothing, I hacked at the script to make it output to STD out before the maize meeting so the -o option is reqired but does not do anything
* This can be modified in the getopts portion of the script, but I left it there in case we wanted a file output
* Takes about 30 minutues to run for one file
```python
python ../QueryDB_CDS.py \
-i 500kb_B73_Mo17_recip_SM.txt  \
-o Arabidopsis_Blast/5500kb_B73_Mo17_recip_SM_Arabidopsis.txt \
-m /home/hirschc1/broha006/software/synmap/b73-phb47/b73.repcds.fa \
-s /panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa \
-d ../../Data_Bases/ArabidopsisRepGenes \
>> Arabidopsis_Blast/500kb_B73_Mo17_recip_SM_Arabidopsis.txt

python ../QueryDB_CDS.py \
-i 500kb_B73_Mo17_recip_SM.txt  \
-o Arabidopsis_Blast/5500kb_B73_Mo17_recip_SM_Arabidopsis.txt \
-m /home/hirschc1/broha006/software/synmap/b73-phb47/b73.repcds.fa \
-s /panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa \
-d ../../Data_Bases/ArabidopsisRepGenes \
>> Arabidopsis_Blast/500kb_B73_Mo17_recip_SM_Arabidopsis.txt
```
* It is also important to mention that I am using the >> to write to file since this script only works in one direction and will complain halfway in the split merged file since gene directions are reversed
* I wrote the commands so both directions would write to the same file then I would anti grep the file for "NOT CALLABLE" to get to the filtered list
```bash
grep -v "NOTCALLABLE" Arabidopsis_Blast/500kb_B73_Mo17_recip_SM_Arabidopsis.txt > Arabidopsis_Blast/500kb_B73_Mo17_recip_SM_Arabidopsis_Filtered.txt

# Here is a for loop bash version for all files in a DIR
for i in $( ls ); do
    NAME=$(basename $i | cut -d. -f1)
    grep -v CALLABLE ${i} > ${NAME}_filtereted.txt
done
```