cd /panfs/roc/data_release/2/umgc/hirschc1/hiseq/180522_D00635_0375_BCCJL5ANXX/Hirsch2_Project_001

#Manually declage the GFF Star indexes as well as the prefixes that we want to use
declare -a INDEX=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_B73_LT"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_BTX_623_LT"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_PH207_LT"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_W22_LT")

declare -a PREFIX=("LT_B73_Ref"
                   "LT_BTX_623_Ref"
                   "LT_PH207_Ref"
                   "LT_W22_Ref")

declare -a GFF=("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Zea_mays.AGPv4.33_Longest_transcript_Exons.gff"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Sbicolor_313_v3.1.gene_exons_Longest_Transcript_Exon.gff3"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/ZmaysPH207_443_v1.1.gene_exons_Longest_Transcript_Exon.gff3"
                  "/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/zea_maysw22_core_fixchr_Longest_Transcript_Exon.gff")

declare -a SCRIPT=("/home/hirschc1/mich0391/Projects/SplitGenes/Scripts/NoCut_Maize_RNA_Counts_LT_Name_B73_W22.sh"
                   "/home/hirschc1/mich0391/Projects/SplitGenes/Scripts/NoCut_Maize_RNA_Counts_LT_ID_PH207_BTX623.sh"
                   "/home/hirschc1/mich0391/Projects/SplitGenes/Scripts/NoCut_Maize_RNA_Counts_LT_ID_PH207_BTX623.sh"
                   "/home/hirschc1/mich0391/Projects/SplitGenes/Scripts/NoCut_Maize_RNA_Counts_LT_Name_B73_W22.sh")

#Grab the headers of each file to pull out the index and output the alignment commands for GNU
for j in {0..3}; do
    for FILE in *.gz; do
        [ -e "$FILE" ] || continue
        OUT="/panfs/roc/scratch/michnoj0/SplitGenes"
        DIR=$PWD
        I7=$(zcat $FILE | head -n1 | cut -f 10 -d":" | cut -f1 -d"+")
        I5=$(zcat $FILE | head -n1 | cut -f 10 -d":" | cut -f1 -d"+")
        #echo "-f $DIR/$file"
        (echo ${SCRIPT[$j]} -f $DIR/$FILE -a $I7 -b $I5 -u ${INDEX[$j]}  -t 24 -g ${GFF[$j]} -p ${PREFIX[$j]} -o $OUT)
    done
done
