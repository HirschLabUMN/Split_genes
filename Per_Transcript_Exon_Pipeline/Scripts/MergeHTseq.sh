# Sctip to merge HTseq files, keep in mind that the line headers will not be present in a file
# and that the order of the $FILES vaiable is the order of your columns. 

mkdir -p /scratch.global/michnoj0/SplitGenes/Counts/Merged
cd /scratch.global/michnoj0/SplitGenes/Counts/Merged

FILES=$(ls ../*B73_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_B73_Ref_HTseq.txt

FILES=$(ls ../*PH207_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_PH207_Ref_HTseq.txt

FILES=$(ls ../*W22_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_W22_Ref_HTseq.txt

FILES=$(ls ../*BTX_623_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_BTX_623_Ref_HTseq.txt

mkdir -p /panfs/roc/scratch/michnoj0/SplitGenes/Counts_LT/Merged
cd /panfs/roc/scratch/michnoj0/SplitGenes/Counts_LT/Merged

FILES=$(ls ../*B73_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_LT_B73_Ref_HTseq.txt

FILES=$(ls ../*PH207_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_LT_PH207_Ref_HTseq.txt

FILES=$(ls ../*W22_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_LT_W22_Ref_HTseq.txt

FILES=$(ls ../*BTX_623_Ref*);
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > Samples_LT_BTX_623_Ref_HTseq.txt
