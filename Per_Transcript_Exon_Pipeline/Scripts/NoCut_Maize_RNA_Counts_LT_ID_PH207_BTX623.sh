#!/bin/bash

#This script takes in your forward reads, reverse reads, 6bp adapter, and output directory and returns
##all intermitant files with a sorted bam being the final output
#This script is to version two of the genome.


set -euo pipefail
#use shebang to indicate the path to the interpretor used to execute this script
#set -e terminated the script if any command exits with a non-zero value
#set -u will prevent bash from running if there is an unset variable
#set -o pipefail This setting prevents errors in a pipeline from being masked. 
#If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.

###### FOR MSI USE
#Please modify the module loads commands to the appropritate software versions
# module load fastqc
# module load bowtie2
# module load bwa
# module load samtools/1.5
# module load bedtools
# module load fastqc
# module load cutadapt
# module load java
# SILENT_JAVA_OPTIONS="$_JAVA_OPTIONS"
# unset _JAVA_OPTIONS

#use getopts for command line arguments
while getopts hf:r::a::b::m::t::u::g:p::o: flag; do
    case $flag in
      # this is the help command
        h)
            echo "This is your help information:

            This script was designed for use on MSI by the Stupar Lab.
            Some of these file paths are hard coded for use on MSI but
            can be easily changed using parameters.

            Usage:
            Maize_RNA_Counts.sh -f forward_read.fastq -r <reverse_read.fastq> -o output_directory 
            
            Version: 1.0

            Options:
                    -a i7 adapter barcode
                    -b i5 adapter barcode
                    -m True :If flag is set to 'True', then hisat2 will
                        be used (default STAR).
                    -t INT :number of threads (default 2).
                    -u Set the location to your STAR index DIR or .fa
                        reference file for HISAT2.
                    -g Path to your genomes GFF file.
                    -p Prefix to end samplename with
                    -o output DIR.
                    "
            exit 2
            ;;
        f)
            echo "Forward read is: $OPTARG";
            FILE=$OPTARG
            I=1
            A="False"
            ALIGNER="STAR"
            THREADS=2
            # hard coded paths are for me, will be overwritten
                # if other users use the corresponding flag.
            REFERENCE="/panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa"
            GFF="/panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/annotation/"
            PREFIX=""
            ;;
        r)
            echo "Reverse read is: $OPTARG";
            FILE2=$OPTARG
            I=2
            ;;
        a)
            echo "Adapter i7 sequence is: $OPTARG";
            ADAPTER=$OPTARG
            A="True"
            ;;
        b)
            echo "Adapter2 i5 sequence is: $OPTARG";
            ADAPTER2=$OPTARG
            A="True"
            ;;
        m)
            echo "You are opting to use HISAT2 instead of STAR";
            ALIGNER="HISAT2"
            ;;
        t)
            echo "$OPTARG threads will be used";
            THREADS=$OPTARG
            ;;
        u)
            echo "$OPTARG is your STAR ref or HISAT Ref";
            REFERENCE=$OPTARG
            ;;
        g)
            echo "$OPTARG is your GFF file";
            GFF=$OPTARG
            ;;
        p)
            echo "Your sample names will ends with $OPTARG ";
            PREFIX=$OPTARG
            ;;
        o)
            echo "everything will be output to this directory: ${OPTARG%/}"
            OUTPUTDIR="${OPTARG%/}"
            ;;
        ?)
            echo "not a valid argument, try using the -h option for usage information"
            exit 2
            ;;
    esac
done

#Letting the user know if the program will use SE or PE pipeline
if [ $I == 2 ]
    then
    echo "setting up for paried-end read pipeline"
    elif [ $I == 1 ]
    then
    echo "setting up for single-end read pipeline"
else
    echo "error, you might be missing a variable try the -h command for help"
    exit 2
fi




#just get the last part of the file name
filename=$(basename "$FILE")
#get what the extension is
extension="${filename##*.}"
#and now you just have the base name for renaming
basename="${filename%%.*}"

#This command will split the base name and grab the first part of it
#It is currently splitting on "_"
#Ex: "Sample1_R1" ==> "Sample1"
samplename=$(echo ${basename} | cut -f1 -d"_" | awk ' {print $1"_"var}' var="$PREFIX" )



#Quality Checks before running the scripts

#check to see that the forward read exists and is readable
if [ -r $FILE ]
    then
    echo "Forward file exists and it readable"
else
    echo "$FILE is not valid, check to see that is a readable fastq file"
    exit
fi

if [ $REFERENCE == "STAR" ]
    then
    if [ -d $REFERENCE ]
        then
        echo "STAR reference DIR is valid"
        else
        echo "This is not a valid directory, check to see that the STAR DIR path is correct"
        exit
    fi
    else
        if [ -r $REFERENCE ]
            then
            echo "Reference .fa file is readable"
        else
            echo "Reference .fa file is not readable"
            exit
        fi
fi



#check same for reverse read if option is flagged
if [ $I == 2 ]
   then
   if [ -r $FILE2 ]
      then
            echo "Forward file exists and it readable"
            filename2=$(basename "$FILE2")
            #get what the extension is
            extension2="${filename2##*.}"
            #and now you just have the base name for renaming
            basename2="${filename2%%.*}"
            #basename2="${filename2%%.anqrp.fastq.gz}"
   else
        echo "$FILE2 is not valid, check to see that is a readable fastq file"
        exit
   fi
fi


if [ -d $OUTPUTDIR  ]; then
   echo "$OUTPUTDIR is a valid directory"
else
   echo "$OUTPUTDIR is not directory"
   exit
fi


#go to Output directory
cd "${OUTPUTDIR}"

#output everything into a log file for easy error identification
echo "Your log file will be located in ${OUTPUTDIR}/FQ2RM_${samplename}.log"
exec > >(tee "${OUTPUTDIR}/FQ2FPKM_${samplename}.log") 2>&1


# #load and run fastqc module for forward and reverse at the same time
# if [ $I == 2 ]
#     then
#         #Run both fastqc programs at once
#         parallel --gnu -j 2 \
#         fastqc -f fastq -o ${OUTPUTDIR} \
#         {1} \
#         ::: ${FILE} ${FILE2}
# else
#     fastqc -f fastq -o ${OUTPUTDIR} ${FILE}
#     echo "Initial FastQC on forward reads worked"
# fi

# unset _JAVA_OPTIONS

# #load and run cutadapt
# echo "Starting adapter trimming"
# REVERSEADAPTER2=$(echo AATGATACGGCGACCACCGAGATCTACAC${ADAPTER2}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG | rev | tr ATGC TACG)
# if [ $A == "True" ]
#     then
#     REVERSEADAPTER1=$(echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}CGATCTCGTATGCCGTCTTCTGCTT" | rev | tr ATGC TACG)
#     if [ $I == 2 ]
#        then
#           #load and run cutadapt for paired end with adaptor
#           cutadapt \
#           -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}CGATCTCGTATGCCGTCTTCTGCTT \
#           -b AATGATACGGCGACCACCGAGATCTACAC${ADAPTER2}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#           -B ${REVERSEADAPTER1} \
#           -B ${REVERSEADAPTER2} \
#           -f fastq \
#           -m 30 \
#           -q 20 \
#           --quality-base=33 \
#           -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
#           -p ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
#           ${FILE} \
#           ${FILE2} \
#           > ${OUTPUTDIR}/${basename}_cutadapt.log

#           #run Fastqc to make sure reads were trimmed
#           fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq &
#           fastqc -f fastq ${OUTPUTDIR}/${basename2}_cutadapt.fastq 

#        else
#           #load and run cutadapt for single end reads with 6bp adaptor
#           cutadapt \
#           -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}CGATCTCGTATGCCGTCTTCTGCTT \
#           -b AATGATACGGCGACCACCGAGATCTACAC${ADAPTER2}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#           -b GCATACGAGATAGCTTCAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT \
#           -f fastq \
#           -m 30 \
#           -q 20 \
#           --quality-base=33 \
#           -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
#           ${FILE} \
#           > ${OUTPUTDIR}/${basename}_cutadapt.log

#           fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq
#     fi
# fi

# #if no adapter is applied then just run cutadadapt with standard illumina
# #to trim based on quiality scores regardless
# if [ $A == "False" ]
#   then
#   echo "Trimming based on quality scores"
#       if [ $I == 2 ]
#         then
#           #load and run cutadapt just searching for standard illumina adapters
#           cutadapt \
#           -b AATGATACGGCGACCACCGAGATCTACAC${ADAPTER2}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#           -B ${REVERSEADAPTER2} \
#           -f fastq \
#           -m 30 \
#           -q 20 \
#           --quality-base=33 \
#           -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
#           -p ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
#           ${FILE} \
#           ${FILE2} \
#           > ${OUTPUTDIR}/${basename}_cutadapt.log

#           #run Fastqc to make sure reads were trimmed
#           fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq &
#           fastqc -f fastq ${OUTPUTDIR}/${basename2}_cutadapt.fastq 
#           wait
#         else
#           #load and run cutadapt just searching for standard illumina adapters
#           cutadapt \
#           -b AATGATACGGCGACCACCGAGATCTACAC${ADAPTER2}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#           -f fastq \
#           -m 30 \
#           -q 20 \
#           --quality-base=33 \
#           -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
#           ${FILE} \
#           > ${OUTPUTDIR}/${basename}_cutadapt.log

#           fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq
#       fi
# fi

# echo "Adapter trimming finished"


#this is the default aligner
if [ $ALIGNER == "STAR" ]
    then
    echo "starting STAR alignment"
    if [ $I == 2 ]
        then
        STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${REFERENCE} \
        --sjdbGTFfile ${GFF} \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon exon \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --readFilesIn ${OUTPUTDIR}/${basename}_cutadapt.fastq ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
        --outFileNamePrefix STAR_${samplename}
        #convert the sam file to a bam file
        samtools view -@ ${THREADS} -bSq 2 ${OUTPUTDIR}/STAR_${samplename}Aligned.out.sam > ${OUTPUTDIR}/STAR_${samplename}.bam
        #we made the bam and no longer need the sam
        rm STAR_${samplename}Aligned.out.sam
        #sort and index the bam file
        samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/STAR_${samplename}.sorted.bam ${OUTPUTDIR}/STAR_${samplename}.bam
        rm ${OUTPUTDIR}/STAR_${samplename}.bam
        samtools index -@ ${THREADS} ${OUTPUTDIR}/STAR_${samplename}.sorted.bam
        echo "STAR complete"
    else
        STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${REFERENCE} \
        --sjdbGTFfile ${GFF} \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon exon \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --readFilesIn ${OUTPUTDIR}/${basename}_cutadapt.fastq \
        --outFileNamePrefix STAR_${samplename}

        #convert the sam file to a bam file
        samtools view -@ ${THREADS} -bSq2 ${OUTPUTDIR}/STAR_${samplename}Aligned.out.sam > ${OUTPUTDIR}/STAR_${samplename}.bam
        #we made the bam and no longer need the sam
        rm STAR_${samplename}Aligned.out.sam
        #sort and index the bam file
        samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/STAR_${samplename}.sorted.bam ${OUTPUTDIR}/STAR_${samplename}.bam
        rm ${OUTPUTDIR}/STAR_${samplename}.bam
        samtools index -@ ${THREADS} ${OUTPUTDIR}/STAR_${samplename}.sorted.bam
        echo "STAR complete"
    fi
fi



if [ $ALIGNER == "HISAT2" ]
    then
    HISATREF=${REFERENCE%.fa}
    echo "starting HISAT alignment"
    if [ $I == 2 ]
        then
        hisat2 \
        -p ${THREADS} \
        --min-intronlen 20 \
        --max-intronlen 20000 \
        -x ${HISATREF} \
        -1 ${OUTPUTDIR}/${basename}_cutadapt.fastq \
        -2 ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
        -S ${OUTPUTDIR}/hisat2_${samplename}.sam \
        #convert the sam file to a bam file
        samtools view -@ ${THREADS} -bSq2 ${OUTPUTDIR}/hisat2_${samplename}.sam > ${OUTPUTDIR}/hisat2_${samplename}.bam
        #we made the bam and no longer need the sam
        rm ${OUTPUTDIR}/hisat2_${samplename}.sam
        #sort and index the bam file
        samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/hisat2_${samplename}.sorted.bam ${OUTPUTDIR}/hisat2_${samplename}.bam
        rm ${OUTPUTDIR}/hisat2_${samplename}.bam
        samtools index -@ ${THREADS} hisat2_${samplename}.sorted.bam
        echo "STAR complete"
    else
        hisat2 \
        -p ${THREADS} \
        --min-intronlen 20 \
        --max-intronlen 20000 \
        -x ${HISATREF} \
        -U ${OUTPUTDIR}/${basename}_cutadapt.fastq \
        -S ${OUTPUTDIR}/hisat2_${samplename}.sam \
        #convert the sam file to a bam file
        samtools view -@ ${THREADS} -bSq2 ${OUTPUTDIR}/hisat2_${samplename}.sam > ${OUTPUTDIR}/hisat2_${samplename}.bam
        #we made the bam and no longer need the sam
        rm ${OUTPUTDIR}/hisat2_${samplename}.sam
        #sort and index the bam file
        samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/hisat2_${samplename}.sorted.bam ${OUTPUTDIR}/hisat2_${samplename}.bam
        rm ${OUTPUTDIR}/hisat2_${samplename}.bam
        samtools index -@ ${THREADS} hisat2_${samplename}.sorted.bam
        echo "STAR complete"
    fi
fi

mkdir -p Counts_LT
#use cufflinks to generate FPKM values
if [ $ALIGNER == "STAR" ]
  then
  htseq-count \
    -f bam \
    -s no \
    -a 0 \
    -t exon \
    -i ID \
    ${OUTPUTDIR}/STAR_${samplename}.sorted.bam \
    ${GFF} \
    > Counts_LT/${samplename}_counts.txt

  else
  htseq-count \
      -f bam \
      -s no \
      -a 0 \
      -t exon \
      -i ID \
      ${OUTPUTDIR}/hisat2_${samplename}.sorted.bam \
          ${GFF} \
      > Counts_LT/${samplename}_counts.txt

fi



