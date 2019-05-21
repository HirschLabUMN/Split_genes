# Split-gene misannotation pipeline

The objective of the code in this repository is to: 1.) generate a key of homologous genes across multiple maize genome annotations, 2.) identify split-gene misannotations, wherein a single gene is incorrectly annotated as two distinct genes, and 3.) use expression data to resolve the correct gene model.  Accurate and consistent genome annotations are necessary for future research in maize, the most well-studied and economically important crop species in the US.

The first two objectives are accomplished by scripts in the subdirectory "Split_Merge_Pipeline," which contains it's own README detailing the use of those scripts, specifically.

For the third objective, RNAseq data must first be mapped and counted according to the pipeline described in the subdirectory "Per_Transcript_Exon_Pipeline."  With the output of the "Split_Merge_Pipeline," we first run scripts/parse_candidates.py program to prepare this output to be input into the scripts/M2f_classify.R program.  To run the parse_candidates.py program, do:

    python parse_candidates.py \
        -i <SplitMerge_output> # Output from Split_Merge_Pipeline.  Reciprically identified split-gene candidates.
        -o <output_directory>
        -n <SampleName1,SampleName2> # comma-separated list (no spaces) that provides the names of each genotype as labelled in the RNAseq counts
        -a <Annotation1,Annotation2> # comma-separated list (no spaces) that provides the two annotation files being compared.  Must be same order as previous option.
        -t <tandem_dup_threshold> # proportion of overlapping BLAST coordines for split-genes with respect to total covered CDS.  Default is 0.1.  Only split-genes with a lower threshold are retained for further analysis
        -e <minimum_exons> # for generating simulated split-genes.  minimum number of exons in a gene in order to generate simulated split-genes from it.  Default is 4
        -S <split_proportion> # Proportion of genes in annotation from which we want to generate simulated split-genes
        -M <merged_proportion> # Proportion of genes in annotation from which we want to generate simulated merged genes

This program will output two files in the specified output directory, one for each of the annotations that were compared.  For example, upon running the following command:

    python scripts/parse_candidates.py -i input/500kb_B73_W22_recip_SM.txt -o ~/Documents/Research/Maize/Split_genes/input -n W,B -a /Users/pmonnahan/Documents/Research/Maize/MaizeSV/misc/zea_maysw22_core_fixchr_Longest_Transcript_Exon.gff,/Users/pmonnahan/Documents/Research/Maize/MaizeSV/misc/Zea_mays.AGPv4.33_Longest_transcript_Exons.gff

two files are output into /Documents/Research/Maize/Split_genes/input: Bsplits_WB_td0.1_m4.txt and Wsplits_WB_td0.1_m4.txt.  "WB" indicates that the W and B are the two annotations being compared. "Bsplits" indicates that this file contains the split-genes in B that correspond to a single gene in W.  "Wsplits" indicates the opposite.  "td0.1" 

With the output from this program, run the M2f_classify.R script to simulate null distributions and make a classification of each split-gene candidate based on a user-defined percentile of null distributions.  Arguments are position and therefore must be in correct order

    Rscript scripts/M2f_classify.R <splits_file> <counts_file> <output_file> <sample_ID> <upper_percentile> <lower_percentile> <min_expression> <read_length> <num_samps>

where:

    splits_file = output containing the split-genes identified from a particular annotation.  Output from parse_candidates.py.  
    counts_file = output from HTseq featureCounts on a per-exon basis
    output_file = name of output file you wish to create
    sample_id = This should uniquely identify all of the HTseq columns that pertain to the sample corresponding to the split-genes being analyzed. E.g. "W" was sampleID for W22 count data, so "W" would identify the different replicates and tissues corresponding to W22.
    upper_percentile = Upper percentile of the simulated split-gene null distribution. Expressed as a proportion.  Observed M2f values above this percentile will be classified that split-genes should remain annotated as distinct genes
    lower_percentile = Lower percentile of the simulated merged-gene null distribution. Expressed as a proportion.  Observed M2f values below this percentile will be classified that split-genes should be merged into a single gene
    min_expression = Minimum normalized expression for a gene to be included in the M2f analysis. Default is 0.01.
    read_length = Length of RNAseq reads
    num_samps = The number of unique samples (columns in HTseq count file) that correspond to the sampleID. E.g. 2 biological reps, each with 10 tissues for a given genotype = 20 for num_samps