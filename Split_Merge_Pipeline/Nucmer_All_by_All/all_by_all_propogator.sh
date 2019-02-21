#Manually declage the GFF Star indexes as well as the prefixes that we want to use
declare -a GFF=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/B73_genes.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/PH207_genes_noV1.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/W22_genes.gff3"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/Mo17_gene.gff"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/PHB47_gene.gff"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/BTX623_genes.gff3")

declare -a DB=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/B73_CDS"
               "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/PH207_CDS"
               "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/W22_CDS"
               "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/Mo17_CDS"
               "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/PHB47_CDS"
               "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/BTX623_CDS")

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
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Data_Bases/phb47.repcds.fa"
                "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize/CDS_results/BTX_623_LongestTS.fasta")

declare -a SPLIT=("d"
                  "a"
                  "b"
                  "a"
                  "v"
                  "G")

COORDS="/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Mummer_Results"

for QUERY in {0..5}; do
    for SUBJECT in {0..5}; do
        # We don't want to compare the same line to itself
        if [ "$QUERY" -ne "$SUBJECT" ]; then
        (echo python All_by_All_Blast.py \
            -q ${GFF[$QUERY]} \
            -s ${GFF[$SUBJECT]} \
            -b ${DB[$SUBJECT]} \
            -l ${CDS[$QUERY]} \
            -o ${LINES[$QUERY]}_${LINES[$SUBJECT]}_AllbyAll_res.txt \
            -n ${COORDS}/${LINES[$QUERY]}_${LINES[$SUBJECT]}_c1000.fil.coords \
            -g ${SPLIT[$SUBJECT]} )
        fi
    done
done

