#!/bin/bash
set -euo pipefail

### DAGchainer synteny pipeline
declare -a GFF=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/B73_genes.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/PH207_genes_noV1.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/W22_genes.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/Mo17_gene.gff"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/PHB47_gene.gff"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/BTX623_genes.gff3")

declare -a LINES=("B73"
                  "PH207"
                  "W22"
                  "Mo17"
                  "PHB47"
                  "BTX623")

declare -a CDS=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/B73_LongestTS.fasta"
                "/home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/W22_LongestTS.fasta"
                "/panfs/roc/groups/13/stuparr/mich0391/software/synmap/mo17-b73/mo17.repcds.fa"
                "/panfs/roc/groups/13/stuparr/mich0391/software/synmap/b73-phb47/phb47.repcds.fa"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/BTX_623_LongestTS.fasta")

for QUERY in {0..5};do
        # preprocess all of the gff's first to subset genes and add gene order
        (python GFF_process.py -i ${GFF[$QUERY]} > ${LINES[$QUERY]}_DAG.gff)
done

for QUERY in {0..5}; do
    for SUBJECT in {0..5}; do
        # We don't want to compare the same line to itself
        if [ "$QUERY" -ne "$SUBJECT" ]; then
        # Blast CDS to CDS
        (echo blastn -query ${CDS[$QUERY]} -subject ${CDS[$SUBJECT]} -evalue 1e-4 -outfmt 6 -out ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS.txt)

        # Convert file to DAG format
        (echo python Blast_2_DAGformat.py -i ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS.txt -o ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS_DAGfmt.txt -q ${LINES[$QUERY]}_DAG.gff -s ${LINES[$SUBJECT]}_DAG.gff)

        # Filter repetative matches
        (echo perl /panfs/roc/groups/13/stuparr/mich0391/Software/DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 5 '<' ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS_DAGfmt.txt '>' ${LINES[$QUERY]}_${LINES[$SUBJECT]}.filtered)

        # Find Syntenic genes
        (echo perl /panfs/roc/groups/13/stuparr/mich0391/Software/DAGCHAINER/run_DAG_chainer.pl -i ${LINES[$QUERY]}_${LINES[$SUBJECT]}.filtered -Z 12 -D 10 -g 1 -A 5)
        fi
    done
done

# for QUERY in {0..5}; do
#     for SUBJECT in {0..5}; do
#         # We don't want to compare the same line to itself
#         if [ "$QUERY" -ne "$SUBJECT" ]; then
#         # Convert file to DAG format
#         (echo python Blast_2_DAGformat.py -i ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS.txt -o ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS_DAGfmt.txt -q ${LINES[$QUERY]}_DAG.gff -s ${LINES[$SUBJECT]}_DAG.gff)
#         fi
#     done
# done

# for QUERY in {0..5}; do
#     for SUBJECT in {0..5}; do
#         # We don't want to compare the same line to itself
#         if [ "$QUERY" -ne "$SUBJECT" ]; then

#         # Filter repetative matches
#         (echo perl /panfs/roc/groups/13/stuparr/mich0391/Software/DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 5 '<' ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_CDS_DAGfmt.txt '>' ${LINES[$QUERY]}_${LINES[$SUBJECT]}.filtered)
#         fi
#     done
# done

# for QUERY in {0..5}; do
#     for SUBJECT in {0..5}; do
#         # We don't want to compare the same line to itself
#         if [ "$QUERY" -ne "$SUBJECT" ]; then
#         # Find Syntenic genes
#         (echo perl /panfs/roc/groups/13/stuparr/mich0391/Software/DAGCHAINER/run_DAG_chainer.pl -i ${LINES[$QUERY]}_${LINES[$SUBJECT]}.filtered -Z 12 -D 10 -g 1 -A 5)
#         fi
#     done
# done

