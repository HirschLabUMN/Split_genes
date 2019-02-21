'''
    File name: make_fakeSplitMerge.py
    Author: Patrick Monnahan
    Date created: 09/01/18
    Python Version: 3.6
    Project: Split Genes
    Upstream of: calcVarRatios.R
    Downstream of: JMM's longest transcript code
    Description: This program generates fake split and merged sets of genes to be used as a null comparison with putative split/merged genes via calcVarRatios.R
'''

import argparse
import pdb
import random
from math import ceil
import operator

def splitGenes(gene_file,):
    split_list = []
    with open(gene_file, 'r') as splits:
        for line in splits:
            split_list.append(line.strip())
    return(split_list)

def getMate(gene, idx):
    trail = gene[8:] # Retrieve trailing part of geneID that consists of just integers
    new_trail = str(int(trail) + idx)
    N = len(trail) - len(new_trail)
    next_gene = gene[:8] + "".join(["0" for k in range(0,N)]) + new_trail # add gene prefix plus new suffix
    return(next_gene)

def getFakes(annotation_path, split_list, split_prop, merge_prop):
    splits = {}
    merged_A = {}
    merged_B = {}
    unchanged = {}
    with open(annotation_path, 'r') as gff:
        for i, line in enumerate(gff):
            if line[0][0] == "#": continue
            line = line.strip().split("\t")
            if line[2] == "gene":
                gene = line[8].strip().split("=")[1].split("_")[0].split(";")[0].split(".")[0]
                if gene not in split_list: # Only make fake split/merges if this gene is not part of a putative split/merge
                    rx = random.uniform(0, 1)
                    if gene not in merged_B: # Prevents this gene + next_gene from being chosen for merging if this gene was already part of a previous merge
                        if rx < split_prop: # make a random split; p is min number of exons
                            splits[gene] = []
                        elif rx > split_prop and rx < (split_prop + merge_prop): # make a random merge with this gene and the one that follows
                            next_gene = getMate(gene, 1)
                            if next_gene not in split_list: # gene has already been checked for, but still need to make sure next_gene is not part of a putative split/merge
                                merged_A[gene] = []
                                merged_B[next_gene] = []
                        else:
                            unchanged[gene] = []
            else: # Each gene of interest has been stored as key in dict, so we will add all corresponding exons as items of these keys
                info = line[8].strip().split(";")  
                gene = [k.split("=")[1] for k in info if "Parent" in k][0].split("_")[0]
                # pdb.set_trace()
                start = int(line[3])
                stop = int(line[4])                       
                exon = [k.split("=")[1] for k in info if "ID" in k or "Name" in k]
                assert len(exon) == 1
                if gene in splits: splits[gene].append([exon[0], line[0], start, stop])
                elif gene in merged_A: merged_A[gene].append([exon[0], line[0], start, stop]) 
                elif gene in merged_B: merged_B[gene].append([exon[0], line[0], start, stop]) 
                elif gene in unchanged: unchanged[gene].append([exon[0], line[0], start, stop])
    return(splits, merged_A, merged_B, unchanged)


def writeFakes(fake_splits, fake_merged_A, fake_merged_B, unchanged, outfile, min_exons):
    with open(outfile + ".txt", 'w') as out:
        for Gene, exon_list in fake_splits.items():
            num_exons = len(exon_list)
            if num_exons >= min_exons:
                exon_list.sort(key=operator.itemgetter(2)) # Sort exons by start position
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3]
                bound = random.randint(1, num_exons - 1)
                for j, exon in enumerate(exon_list):
                    merged = exon.copy() # This is the original parent
                    # bound = ceil(num_exons / 2)
                    if j < bound: # assign half of the exons to split_a and half to split_b
                        exon += [f"{Gene}a", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "1", "fakeSplit"]
                    else:
                        exon += [f"{Gene}b", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "1", "fakeSplit"]
                    merged += [f"{Gene}", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "1", "fakeSplit"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
                    out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in fake_merged_A.items(): 
            for j, exon in enumerate(exon_list):
                Partner = getMate(Gene, 1)
                merged = exon.copy()
                exon += [Gene, gene_start, gene_end, f"{Gene}c", f"{Gene},{Partner}", "1", "fakeMerged"]
                merged += [f"{Gene}c", gene_start, gene_end, f"{Gene}c", f"{Gene},{Partner}", "1", "fakeMerged"]
                out.write("\t".join([str(x) for x in exon]) + "\n")
                out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in fake_merged_B.items():
            for j, exon in enumerate(exon_list):
                merged = exon.copy()
                Partner = getMate(Gene, -1)
                exon += [Gene, gene_start, gene_end, f"{Partner}c", f"{Partner},{Gene}", "1", "fakeMerged"]
                merged += [f"{Partner}c", gene_start, gene_end, f"{Partner}c", f"{Partner},{Gene}", "1", "fakeMerged"]
                out.write("\t".join([str(x) for x in exon]) + "\n")
                out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in unchanged.items():
            for j, exon in enumerate(exon_list):
                exon += [Gene, gene_start, gene_end, Gene, f"{Gene},{Gene}", "1", "unchanged"]
                out.write("\t".join([str(x) for x in exon]) + "\n")
    return()
                                   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "This program generates fake split and merged sets of genes to be used as a null comparison with putative split/merged genes via calcVarRatios.R")
    parser.add_argument('-a', type=str, metavar='annotation_list', required=True, help='Comma separated string (no spaces) containing paths to annotation files.  annotation files should exons corresponding to a just a single transcript.')
    parser.add_argument('-q', type=str, metavar='query_genes', required=True, help='File with geneID (one per line) for which you wish to retrieve first and last exon')
    parser.add_argument('-o', type=str, metavar='out_suffix_list', required=True, help="comma separated list of output prefixes in SAME order as annotation list")
    parser.add_argument('-m', type=int, metavar='minimum_exons', required=True, help="minimum number of exons required to make a fake split gene")
    parser.add_argument('-S', type=float, metavar='split_proportion', required=False, default=0.2, help="proportion of genes for which we want to make fake splits")
    parser.add_argument('-M', type=float, metavar='merged_proportion', required=False, default=0.3, help="proportion of genes for which we want to make fake merge")
    # parser.add_argument('-s', type=float, metavar='syntenic_genes', required=True, help="File containing syntenic genes for B and P. One gene name per line")
    parser.add_argument('-v', action="store_true")
    args = parser.parse_args()

    gffs = args.a.split(",")
    outs = args.o.split(",")

    split_list = splitGenes(args.q)

    for k, path in enumerate(gffs):
        fake_splits, fake_merged_A, fake_merged_B, unchanged = getFakes(path, split_list, args.S, args.M)
        writeFakes(fake_splits, fake_merged_A, fake_merged_B, unchanged, outs[k], args.m)


