# Split gene per exon bioinformatics pipeline
This README is designed for one of the hirsch lab members to reproduces the 
split gene counts per gene/exon analysis

## Project summary
The purpose of this pipeline is to use RNA-seq from B73, PH207, and W22 and 
quantify expression levels on a per gene and per exon basis to see whiwch genes 
across references are truely split and both expressed, if only one of the two 
are expressed, and only half the combination of the two genes are expressed.

The main DIR where all of the code and logs are stored on MSI:
```
/panfs/roc/groups/14/hirschc1/shared/projects/splits/Per_Transcript_Exon_Pipeline
```

* Any file that has an LT in the name is part of the per exon pipeline for
the longest transcript.


Log DIRs for trimmed reads and mapping statistics
```
Per_Transcript_Exon_Pipeline/Logs/Adapter_Trimming_logs
Per_Transcript_Exon_Pipeline/Logs/Alignment_logs
```

HTSeq individual counts and merged files in to matricies
```
Per_Transcript_Exon_Pipeline/Counts/Individual
Per_Transcript_Exon_Pipeline/Counts/Matrix
```

## Pipeline description
Scripts used for pipeline and job submission can be found at:
```
Per_Transcript_Exon_Pipeline/Scripts
```

Software versions
HTseq and Cutadapt we both installed using a Conda env
```
fastqc 0.11.7
cutadapt 1.16 (not a default module on MSI, I used conda)
STAR 020201 (Not on MSI)
samtools 1.6
htseq 0.10.0 (Not on MSI)
```

In brief:

* Cutdapt ran with the quality set to 20 and minimum length set to 30 bp
* STAR was run using two-pass mode, and outputted SAMstrand field in case we 
wanted to use stringtie
* Samtools was used to convert SAM to BAM while sorting and indexing. Reads 
with a mapping quality of 2 or lower were also filtered out at this step (note 
that his is not reflected in the STAR mapping statistics log)
* htseq-count was run with the stranded option set to no and min quality set to 
0 (since we already pre-filtered).
* The main job script I set up/use is called "Maize_RNA_Counts.sh" if you want 
all of the details of the pipeline.

## Order of scripts to run pipeline (please read section below before running anything)
1. Longest\_Transcript\_GFF\_Generator.sh
3. STAR_Index_job.sh (MSI job)
3. AlignPropogator.sh
4. AlignPropogator_LT.sh
5. BTX_Align.sh (MSI job)
6. Othe\r_3\_Align.sh (MSI job)
7. Longest\_Transcript\_Align.sh (MSI job)
8. MergeHTseq.sh

### What do the scripts actually do

* Longest_Transcript\_GFF\_Generator.sh
    * Dependencies (Parse\_GFF\_LongestTranscript.py)
    * You will need to modify the -o option to point to whatever DIR you want
    * This is a bash script that uses a custom python script 
    (Parse_GFF_LongestTranscript.py)that will go into a GFF file and return
    exon field belonging to the longest transcript per gene. Run the python 
    script directly with the -h option to get more information on how to run it
    to store your data

* STAR\_Index\_job.sh
    * Dependencies (StarIndex.sh)
    * You will need to modify StarIndex.sh to point DIR paths to locations you 
    want and you have to point STAR index commands to your own gff and fastaa
    files
    * This is a MSI job to run all of the STAR indexing. STAR needs about 10 GB
    for every gigabase of genome and I found that the upper limit of the largest
     fasta was about 52 GB.
    * B73 indexing using the exons for the longest transcript needs to have the 
     exon parent transcript defined since the defualt parameter in STAR does not
      coinside with what is in the GFF

* AlignPropogator.sh and AlignPropogator_LT.sh
    * This bash script will print out all the commands you need for the all by 
    all comparison.
    * You will need to modify the bash script to point the appropriate variables 
    to the right items. Write the output of the command to a .txt file for use 
    later on.

* BTX_Align.sh
    * Dependencies (BTX_ReadMapCommands.txt, Maize_RNA_Counts.sh)
    * You will need to make BTX_ReadMapCommands.txt by just running a grep 
    command from the .txt file you created from AlignPropogator.sh. You will 
    also need to point GNU parallel to the BTX_ReadMapCommands.txt in the bash 
    script
    * Something like grep "BTX" ReadMapCommands.txt > BTX_ReadMapCommands.txt
    * WHY DO WE DO THIS???????
        * The scripts to run the actual processing commands in 
        Maize_RNA_Counts.sh uses cut adapt and we don't need to run cutadpt on
        the samples more that once, we can just skip to the trimmed reads and 
        map to other reference genomes. That is why there is Maize_RNA_Counts.sh 
        as well as NoCut_Maize_RNA_Counts.sh. The NoCut version skips the 
        adapter trimming. 
        * The Sorghum genome is also a lot smaller so we can load two STAR 
        genomes on a single MSI node. This is why we use --jobs set to 2 in GNU 
        parallel

* Othe\r_3\_Align.sh (Hard extra step)
    * Dependencies (Ref_to_Ref_ReadMapCommands.txt, NoCut_Maize_RNA_Counts.sh)
    * Simlar to the above script you need to grep from the ReadMapCommands.txt 
    to get everything but the BTX alignments
    * Something like grep -v "BTX" ReadMapCommands.txt > Ref_to_Ref_ReadMapCommands.txt
    * You will need to modify Othe\r_3\_Align.sh to point at the location of 
    your Ref_to_Ref_ReadMapCommands.txt

* Longest\_Transcript\_Align.sh 
    * Dependencies (ReadMapCommands_LT.txt, NoCut_Maize_RNA_Counts_LT_ID_PH207_BTX623.sh, NoCut_Maize_RNA_Counts_LT_Name_B73_W22.sh)
    * Similar to the above scripts except the .txt file is just the output from 
    AlignPropogator_LT.sh. No need to run any grep commands since all of the 
    original read files are trimmed. 
    * Need to modify the bash script to point at the location of your txt file
    * There are two different alignment scripts because the Exon ID nomenclature
    is not the same between reference GFF files.

* MergeHTseq.sh
    * You will need to modify the mkdir and cd commands to the directory where 
    your Counts are located
    * This is just a series of commands to merge single HTseq files into a 
    Matrix
    * IMPORTANT: This will not give you column names in the file. This column 
    names will has to be added after. The order of names would be the same order 
    as if you ran a "ls" command. 


