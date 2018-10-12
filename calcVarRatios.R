

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

require(data.table)
require(dplyr)
require(tidyr)
require(parallel)
require(stringr)
require(lme4)
require(inline)
require(tibble)

outfile = args[4]
outfile2 = args[5]
counts = read.table(args[1], head=T)
splits = read.table(args[2])
fsplits = read.table(args[3])

oldNames = c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10") 
newNames = c("exon","chrom","pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum")

setnames(splits, old = oldNames, new = newNames)
setnames(fsplits, old = c(oldNames,"V11"), new = c(newNames, "source"))
setnames(counts, old = c("Genes"), new = c("exon"))

splits['source']="real"

error.df <- data.frame("Df" = c(-9, -9), "AIC" = c(-9, -9), "BIC" = c(-9, -9), "logLik" = c(-9, -9), "deviance" = c(-9, -9), "Chisq" = c(-9, -9), "pval" = c(-9, -9), "varRatio" = c(-9, -9), "prog" = c(-9, -9), "parent" = c(-9, -9), "status" = c("error", "error"))

getStats <- function(df, Theta){
  full = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene), data=df, family=MASS::negative.binomial(theta=Theta))
  reduced = glmer(value ~ (1|tissue) + (1|sample), data=df, family=MASS::negative.binomial(theta=Theta))
  cc = as.data.frame(VarCorr(full))
  ratio = cc[cc$grp=="gene",4] / sum(cc[,4])
  dd = as.data.frame(anova(full, reduced))[,-7]
  colnames(dd)[7] = "pval"
  dd['varRatio'] = ratio
  dd['prog'] = df$prog[1]
  dd['parent'] = df$parent[1]
  dd['status'] = "Good"
  dd['source'] = df$source[1]
  return(dd)
}

runNB <- function(i){
  Prog = as.character(unique(pp$prog)[i])
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  results = tryCatch({
    getStats(dat, Theta)
  }, warning = function(w) {

    rr = getStats(dat, Theta)
    rr['status'] = "warning"
    return(rr)
  }, error = function(e) {
    errD = cbind(error.df)
    errD$parent = dat$parent[1]
    errD$prog = dat$prog[1]
    errD$source = dat$source[1]
    return(errD)}) 
  return(results)
}

runLMM <- function(i){
  Prog = as.character(unique(pp$prog)[i])
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  results = tryCatch({
    getStatsLMM(dat)
  }, warning = function(w) {
    print(w)
    rr = getStatsLMM(dat)
    rr['status'] = "warning"
    return(rr)
  }, error = function(e) {
    errD = cbind(error.df)
    errD$parent = dat$parent[1]
    errD$prog = dat$prog[1]
    errD$source = dat$source[1]
    return(errD)}) 
  return(results)
}

getStatsLMM <- function(df, Theta){
  full = lmer(log(value+1) ~ (1|tissue) + (1|sample) + (1 | gene), data=df)
  reduced = lmer(log(value+1) ~ (1|tissue) + (1|sample), data=df)
  cc = as.data.frame(VarCorr(full))
  ratio = cc[cc$grp=="gene",4] / sum(cc[,4])
  dd = as.data.frame(anova(full, reduced))[,-7]
  colnames(dd)[7] = "pval"
  dd['varRatio'] = ratio
  dd['prog'] = df$prog[1]
  dd['parent'] = df$parent[1]
  dd['status'] = "Good"
  dd['source'] = df$source[1]
  return(dd)
}

runOne = function(i){
  Prog = as.character(unique(pp$prog)[i])
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  aa = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene), data=dat, family=MASS::negative.binomial(theta=Theta))
  return(aa)
}

munge = function(counts, realSplits, fakeSplits){
  splits = rbind(realSplits, fakeSplits)
  d = merge(splits, counts, by = c("exon"), all.y=TRUE)
  d.m = melt(d, id.vars = c("exon","chrom" ,"pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum", "source"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  d.n$tissue = as.factor(d.n$tissue)
  d.n$sample = as.factor(d.n$sample)
  d.n$rep = as.factor(d.n$rep)
  d.n = d.n %>% filter(!exon %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual", "__alignment_not_unique")) %>% filter(nlevels(droplevels(gene)) > 1) %>% as.data.frame()
  return(d.n)
}

munge2 = function(fpkm){
  d.m = melt(fpkm, id.vars = c("gene","parent","prog", "famnum", "source","pos.gene","end.gene"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  d.n$tissue = as.factor(d.n$tissue)
  d.n$sample = as.factor(d.n$sample)
  d.n$rep = as.factor(d.n$rep)
  d.n$gene = as.factor(d.n$gene)
  d.n = d.n %>% filter(nlevels(droplevels(gene)) > 1 & source != "unchanged") %>% as.data.frame()
  return(d.n)
}

munge3 = function(tpm_exons){
  d.m = melt(tpm_exons, id.vars = c("exon","gene","parent","prog", "famnum", "source","pos.exon","end.exon","rowid","basepairs","pos.gene","end.gene"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  d.n$tissue = as.factor(d.n$tissue)
  d.n$sample = as.factor(d.n$sample)
  d.n$rep = as.factor(d.n$rep)
  d.n$gene = as.factor(d.n$gene)
  d.n = d.n %>% filter(nlevels(droplevels(gene)) > 1 & source != "unchanged") %>% as.data.frame()
  return(d.n)
}

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

print("Munging...")
p = munge(counts, splits, fsplits)
print("Munging...done")

p['basepairs'] = abs(p$end.exon - p$pos.exon)
gene_lengths = p %>% filter(sample=="B" & tissue=="A" & rep=="R1") %>% group_by(gene) %>% summarize(basepairs = sum(basepairs)) %>% as.data.frame()
#p['stdval'] = p$value / p$basepairs
#p = p %>% group_by(prog, gene) %>% filter(stdval != max(stdval) | stdval != min(stdval)) %>% as.data.frame()

nn = colnames(counts)[2:61]
mm = c(rep("B",20), rep("P",20), rep("W",20))
tt = c(rep(c("A","A","Em","Em","En","En","IE","IE","I","I","L10","L10","L","L","R","R","SC","SC","T","T"),3))
rr = c(rep(c("R1","R2"),30))
coldata = data.frame(row.names = nn, geno = mm, tissue = tt, rep = rr)

# print("Casting count data...")
# nsg = dcast(p, gene + parent + prog ~ variable, value.var = "value", fun.aggregate = sum)
# nsg = nsg %>% filter(!is.na(gene) & !duplicated(gene)) %>% as.data.frame()
# print("Casting count data...done")
# 
# nsg = merge(nsg, gene_lengths, by="gene")
# row.names(nsg) = nsg$gene
# print(head(nsg))
# dds = DESeqDataSetFromMatrix(countData = nsg[,3:62], colData = coldata, design = ~ geno + tissue + rep)
# rownames(dds) = nsg$gene
#
# mcols(dds) <- cbind(mcols(dds), nsg[63])
# dds = DESeq(dds)

# tpm = as.data.frame(counts_to_tpm(nsg[,4:63],nsg$basepairs, c(rep(50,60))))

nse = dcast(p, exon + gene + parent + prog + basepairs + pos.gene + end.gene + pos.exon + end.exon + famnum + source ~ variable, value.var = "value")
nse = nse[nse$basepairs>50,]
nse=nse[!is.na(nse$basepairs),]
nse['rowid'] = paste(nse$exon, nse$prog, nse$gene, nse$parent)
row.names(nse) = nse$rowid
tpm_e = as.data.frame(counts_to_tpm(nse[,12:71], nse$basepairs, c(rep(50,60))))
print(nse)
# qq = p %>% group_by(gene, parent, prog) %>% summarize(pos.gene = min(pos.exon), end.gene = max(end.exon), famnum=min(famnum), source = min(source)) %>% as.data.frame()

# tpm = rownames_to_column(tpm, var = "gene")
tpm_e = rownames_to_column(tpm_e, var = "rowid")

print(head(tpm_e))
# ww = merge(tpm, qq, by=c("gene"), all.x=T)
ww = merge(tpm_e, nse[,c("exon","gene","parent","prog","basepairs","pos.gene","end.gene","pos.exon","end.exon","famnum","source","rowid")], by="rowid", all.x=T)
# print("Extracting DESeq2 data...done")
print(head(ww))
print("Re-munging...")
# pp = munge2(ww)
pp = munge3(ww)
print("Re-munging...done")
print(head(pp))

# Filters?
#   exons minimum tpm 
#   expressed in minimum number of tissues
#   
#   plot varRatio versus mean tpm
#   filter genes with huge CV??
#   
#   
pp = pp %>% group_by(exon) %>% filter(mean(value)>0.01) %>% as.data.frame()
pp = pp %>% group_by(parent) %>% filter(end.exon==max(end.exon) | pos.exon==min(pos.exon)) %>% as.data.frame()

# throw out first and last exon??
# print("Merging dispersion estimates...")
# pp = merge(ee, disp, by=c("gene"), all.x = TRUE)
# print("Merging dispersion estimates...done")

print("Filtering for valid progeny sets...")
# pp = pp %>% filter(disp != 60) %>% as.data.frame() ## FILTER FOR GENES WITH DISPERSION ==60; MEANS VERY LOW COUNTS
# pp = pp %>% group_by(parent, prog) %>% filter(nlevels(droplevels(gene)) > 2) %>% ungroup() %>% as.data.frame()
# print("Filtering for valid progeny sets...done")

print("Writing input data...")
write.table(pp, outfile2, row.names = F, quote=F)
print("Writing input data...done")

# rm(ee)
rm(nse)
rm(ww)
# rm(dds)
# rm(qq)
rm(p)

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

results = data.frame("Df" = as.numeric(), "AIC" = as.numeric(), "BIC" = as.numeric(), "logLik" = as.numeric(), "deviance" = as.numeric(), "Chisq" = as.numeric(), "Chi Df" = as.numeric(), "Pr(>Chisq)" = as.numeric(), "varRatio" = as.numeric(), "prog" = as.character(), "parent" = as.character(), "status" = as.character(), "source" = as.character())

for (i in 1:(length(unique(pp$prog))/10)){
  # for (i in 1:10){
  j = i * 10
  progress = round(i/(length(unique(pp$prog))/10), digits = 4) * 100
  print(paste0("Running jobs: ", j - 9, "-", j, "; (", progress, "% complete)"))
  ff = mclapply((j - 9):j, FUN = runLMM, mc.cores=2)
  wait()
  ff = as.data.frame(do.call(rbind,ff))
  results = rbind(results, ff)
}

write.table(results, outfile, row.names = F, quote=F)
