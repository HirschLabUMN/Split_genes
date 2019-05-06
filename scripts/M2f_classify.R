#################### BEGIN: Load Packages ###############
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(stringr)
library(tibble)
library(magrittr)
library(wrapr)

#mask these functions from other packages
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarize <- dplyr::summarize

##################### END: Load packages #####################


##################### BEGIN: Read arguments ##################

args = commandArgs(trailingOnly=TRUE)

split_file = args[1] # Output from parse_candidates.py
count_file = args[2] # Output from HTseq counts (per exon)
outfile = args[3] 
sampleID = args[4] # String that identifies samples to be used for M2f .e.g. "W" for W22 in our data
upper_quantile = as.numeric(args[5]) # M2f quantile of null distribution above which candidates will be classified as split
lower_quantile = as.numeric(args[6]) # M2f quantile of null distribution below which candidates will be classified as merged
minTPM = as.numeric(args[7]) # Only include a gene in M2f calculation if normalized expression (TPM) is above this value.  
if(minTPM == -9){ # default value of minTPM
  minTPM = 0.01
}
read_length = as.numeric(args[8]) # Needed for normalization.  Length of RNAseq reads
num_obs = as.numeric(args[9]) # Number of observations for the sample in question.  E,g. Number of tissues * number of reps
##################### END: Read arguments ##################

##################### BEGIN: Load Data ##################
splits = read.table(split_file, head=T)
counts = read.table(count_file, head=T)

#Remove any duplicates
splits %<>% distinct()
##################### END: Load Data ##################


########################## BEGIN: Define functions #####################################
# Convert count data to transcripts per million (tpm)
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

# Calculate the average log fold change across a set of potential split genes -- called by calcM2f
getAvgLog = function(v, verbose=F){
  vm = matrix(v, nrow = length(v), ncol = length(v)) # creates a square matrix out of a vector of tpm values
  svm = sweep(vm, 2, v, FUN = "/") # divide rows of matrix by vector (v)...basically get fold change of all relative to all
  lt = abs(log(svm[lower.tri(svm)], 2)) # get upper and lower triangular matrix...i.e. exclude the diagonals of matrix which is just all 1's
  ut = abs(log(svm[upper.tri(svm)], 2)) # abs(log(1/5)) == abs(log(5)).  each fold change gets double counted but this comes out in the average
  m2f = mean(c(lt,ut)) 
  if (verbose==TRUE){
    print(vm)
    print(svm)
    print(lt)
    print(ut)
    print(m2f)
  }
  return(m2f)
}

convertFactors = function(df){
  df$tissue = as.factor(df$tissue)
  df$sample = as.factor(df$sample)
  df$rep = as.factor(df$rep)
  df$gene = as.factor(df$gene)
  return(df)
}

# Parse the sample names in columns of counts file
parseName = function(df){
  df = df %>%mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>%mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  return(df)
}

# Melts combined dataset and labels rows with information in sample_ids
munge = function(df, melt_by){
  d.m = melt(df, id.vars = melt_by)
  d.n = parseName(d.m)
  d.n = convertFactors(d.n)
  return(d.n)
}

#Major function that takes exon-based count data from HTseq (counts), and a file from parse_candidates.py.  If evaluating split genes in W, use count data from W mapped to W
formatData = function(counts, splits, minTPM = 0, sampleID = c("B", "P", "W"), minGenes = 1, fmt="de", readLength = read_length, numObs = num_obs){
  print("Munging...")
  splits %<>% distinct()
  d = merge(splits, counts, by = c("exon"), all.y=TRUE)
  Names = c("exon","pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog", "source")
  p = munge(d, c(Names, "chrom"))
  #Remove entries without unique exon assignment and retain only the samples that we specify
  p %<>% filter(!exon %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual", "__alignment_not_unique")) %>% as.data.frame()
  print("Munging...done")
  p['basepairs'] = abs(p$end.exon - p$pos.exon)
  if (fmt == "m2f"){
    p %<>% filter(sample %in% sampleID) %>% as.data.frame()
    print("Calculating TPM...")
    nse = dcast(p, exon + gene + parent + prog + basepairs + pos.gene + end.gene + pos.exon + end.exon + source ~ variable, value.var = "value")
    nse = nse[nse$basepairs > read_length,]
    nse=nse[!is.na(nse$basepairs),]
    nse['rowid'] = paste(nse$exon, nse$prog, nse$gene, nse$parent)
    row.names(nse) = nse$rowid
    nse %<>% select(matches(paste(paste(c(Names, "basepairs", "rowid"), collapse = "|"), paste(sampleID, ".", collapse = "|", sep = ""), sep = "|")))
    #20 is tissue number * rep number, so total number of counts per sample
    tpm_e = as.data.frame(counts_to_tpm(nse[, 12:ncol(nse)-1], nse$basepairs, c(rep(readLength, num_obs * length(sampleID)))))
    tpm_e = rownames_to_column(tpm_e, var = "rowid")
    ww = merge(tpm_e, nse[,c(Names, "basepairs", "rowid")], by="rowid", all.x=T)
    print("Done calculating TPM...")
    print("Re-munging...")
    p = munge(ww, c(Names, "basepairs", "rowid"))
    #Average TPM across exons and label parents
    p = p %>% group_by(gene, tissue, sample, rep, source, parent, prog, pos.gene, end.gene) %>% summarize(value = mean(value)) %>% as.data.frame()
    p %<>% mutate(isParent = ifelse(as.character(gene) == as.character(parent), TRUE, FALSE))
  }
  return(p)
}

# Calculates an M2f value for all 
calcM2f = function(df, sampleID, minTPM = 0, FUN=mean){
  # Only consider the split/child genes, filter out 0s, unwanted samples
  DF = df %>%filter(isParent == FALSE & value > minTPM & sample == sampleID) %>% group_by(parent, prog, source, sample, rep, tissue) %>% summarize(m2f = getAvgLog(value)) %>% ungroup() %>% group_by(parent, prog, source) %>% summarize(M2f = FUN(m2f,na.rm=T), numTiss=length(unique(as.character(tissue))), BiDiffMean = NaN)
  badPars = DF %>% filter(M2f=="NaN")
  DF2 = df %>%filter(parent %in% badPars$parent & isParent == FALSE & sample == sampleID) %>% mutate(expressed = ifelse(value > 0, 1, 0)) %>% group_by(parent, prog, source, sample, rep, tissue) %>% summarize(diff = mean(abs(diff(expressed)))) %>% ungroup() %>% group_by(parent, prog, source) %>% summarize(M2f = NaN, numTiss=length(unique(as.character(tissue))), BiDiffMean = sum(diff))
  DF %<>% filter(M2f!=0) %>% as.data.frame()
  DF = rbind(DF, as.data.frame(DF2))
  return(DF)
}

#Adds a column to the data frame labelling whether each call in calcM2f results exceed "significance threshold"
labelExceeds = function(df, merged_quantile, split_quantile, Source = c("real")){
  mq = quantile(df[df$source == "simMerged",]$M2f, merged_quantile, na.rm = T)
  sq = quantile(df[df$source == "simSplit",]$M2f, split_quantile, na.rm = T)
  df %<>%filter(source %in% Source) %>%mutate(Output_call = ifelse(M2f > sq, "SPLIT", ifelse(M2f < mq, "MERGED", "NOCALL"))) %>% as.data.frame() 
  return(df)
}
######################## END: DEFINE FUNCTIONS ######################################

######################## BEGIN: CALC M2F ############################################
# HOW TO HANDLE THE SAMPLEID'S IN AN AUTOMATED WAY?
dat = formatData(counts, splits, sampleID = sampleID, fmt = "m2f")

# Calculate the mean log-2-fold expression across splitGene candidates for each set of candidates
M2f = calcM2f(dat, sampleID, minTPM = minTPM)

split_cutoff = quantile(M2f[M2f$source=="simSplit",]$M2f, upper_quantile, na.rm=T)
merge_cutoff = quantile(M2f[M2f$source == "simMerged",]$M2f, lower_quantile, na.rm=T)

M2f %<>% mutate(Call = ifelse(M2f > split_cutoff, "Split", ifelse(M2f < merge_cutoff, "Merged", "NoCall")), ref = sampleID)

write.table(M2f, outfile, quote = F, row.names = F)
######################## END: CALC M2F ############################################

