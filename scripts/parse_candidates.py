'''
    File name: parse_candidates.py
    Author: Patrick Monnahan
    Date created: 04/29/19
    Python Version: 3.6
    Project: M2f
    Upstream of: M2f_classify.R
    Downstream of: Syntenic homology pipeline
    Description: This program generates simulated split and merged sets of genes to be used as a null comparison with putative split/merged genes as identified in the syntenic homology pipeline.  This produces the formatted file that is input alongside the expression information into the M2f_classify.R script
'''

import argparse
import pdb
import random
from math import ceil
import operator
import itertools
import os

def parseFile(in_file, td_thresh):

    def getOverlap(a, b):
        a.sort()
        b.sort()
        return max(0, min(int(a[1]), int(b[1])) - max(int(a[0]), int(b[0]))) / max(0, max(int(a[1]), int(b[1])) - min(int(a[0]), int(b[0])))

    split_dict = {}
    sets = 0
    Ovlps = 0
    NOvlps = 0
    AOP = 0 # average proportion of overlap
    os = 0 # counter for opposite strands
    nos = 0 
    tds = 0 
    nadj = 0
    with open(in_file, 'r') as inFile:
        for i, line in enumerate(inFile):
            if i > 0 and line[0] != "#" and "combinations" not in line:
                line = line.split()
                td = float(line[-1])

                opposite_strands = line[17]
                # print(td, opposite_strands)
                if "non-adjacent_syntenic" not in line and td <= td_thresh and opposite_strands == "False":
                    if opposite_strands == "False": nos += 1
                    sets += 1
                    splits = line[14].split(";")
                    prog_string = ""
                    Coords = []
                    for split in splits:
                        coords = split.split(",")
                        coords = [x.replace('"','') for x in coords]
                        prog_string += coords[0] + ","
                        Coords.append(coords[2:4])
                    ovlps = []
                    for pair in itertools.combinations(Coords, 2):
                        # pdb.set_trace()
                        if getOverlap(*pair) > 0:
                            ovlps.append(pair)
                            # print(splits)
                            Ovlps += 1
                            AOP += getOverlap(*pair)
                        else: NOvlps += 1
                    ovlps = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(ovlps))))
                    #Get blacklisted index of genes that overlap...or just use coordinates and hope they are unique??
                    for split in splits:
                        coords = split.split(",")
                        if coords[2] not in ovlps and coords[3] not in ovlps:
                            out = [coords[1], coords[2], coords[3], coords[0], line[1], prog_string[:-1]]
                            split_dict[coords[0]] = out
                elif "non-adjacent_syntenic" in line: nadj += 1
                elif opposite_strands == "True": os += 1
                elif td >= td_thresh: tds += 1
                
                        # else: pdb.set_trace()
    if Ovlps > 0:
        print(f"# overlaps = {Ovlps}; # non-overlaps = {NOvlps}; prop overlaps {Ovlps / (Ovlps + NOvlps)}; avg. overlap {AOP / Ovlps}")
    else: print(f"no overlaps")
    print(f"# opp. strands = {os}; # non-opp. strands = {nos}; prop opp. strands {os / (os + nos)}")
    print(f"# tandem dups = {tds}; # non-adjacent_syntenic genes = {nadj}")
    print(f"{i} total one-to-many homologies")

    return(split_dict)

def getMate(gene, idx):
    """Used when simulating merged genes.  A gene has been randomly selected and this function gets the next proximal gene"""
    trail = gene[8:] # Retrieve trailing part of geneID that consists of just integers
    new_trail = str(int(trail) + idx)
    N = len(trail) - len(new_trail)
    next_gene = gene[:8] + "".join(["0" for k in range(0,N)]) + new_trail # add gene prefix plus new suffix
    return(next_gene)

def getExons(annotation_path, real_dict, split_prop, merge_prop):
    splits = {}
    merged_A = {}
    merged_B = {}
    unchanged = {}
    Real = {}

    #Get candidate split genes to screen out of simulations
    ss = list(real_dict.keys())
    for s in real_dict.keys():
        ss.append(real_dict[s])

    with open(annotation_path, 'r') as gff:
        for i, line in enumerate(gff):
            if line[0][0] == "#": continue
            line = line.strip().split("\t")
            if line[2] == "gene":
                gene = line[8].strip().split("=")[1].split("_")[0].split(";")[0].split(".")[0] #This is not very generalized
                if gene not in ss: # Only make fake split/merges if this gene is not part of a putative split/merge
                    rx = random.uniform(0, 1)
                    if gene not in merged_B: # Prevents this gene + next_gene from being chosen for merging if this gene was already part of a previous merge
                        if rx < split_prop: # make a random split; p is min number of exons
                            splits[gene] = []
                        elif rx > split_prop and rx < (split_prop + merge_prop): # make a random merge with this gene and the one that follows
                            next_gene = getMate(gene, 1)
                            if next_gene not in real_dict.keys(): # gene has already been checked for, but still need to make sure next_gene is not part of a putative split/merge
                                merged_A[gene] = []
                                merged_B[next_gene] = []
                        else:
                            unchanged[gene] = []
                else:
                    Real[gene] = []
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
                elif gene in Real: Real[gene].append([exon[0], line[0], start, stop])
    return(splits, merged_A, merged_B, unchanged, Real)

def writeSim(sim_splits, sim_merged_A, sim_merged_B, unchanged, real, real_info, outfile, min_exons):
    with open(outfile + ".txt", 'w') as out:
        out.write("exon\tchrom\tpos.exon\tend.exon\tgene\tpos.gene\tend.gene\tparent\tprog\tsource\n")
        for Gene, exon_list in sim_splits.items():
            num_exons = len(exon_list)
            if num_exons >= min_exons:
                exon_list.sort(key=operator.itemgetter(2)) # Sort exons by start position
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3]
                bound = random.randint(1, num_exons - 1)
                for j, exon in enumerate(exon_list):
                    merged = exon.copy() # This is the original parent
                    if j < bound: # assign half of the exons to split_a and half to split_b
                        exon += [f"{Gene}a", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "simSplit"]
                    else:
                        exon += [f"{Gene}b", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "simSplit"]
                    merged += [f"{Gene}", gene_start, gene_end, Gene, f"{Gene}a,{Gene}b", "simSplit"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
                    out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in sim_merged_A.items():
            if len(exon_list)>0:
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3] 
                for j, exon in enumerate(exon_list):
                    Partner = getMate(Gene, 1)
                    merged = exon.copy()
                    exon += [Gene, gene_start, gene_end, f"{Gene}c", f"{Gene},{Partner}", "simMerged"]
                    merged += [f"{Gene}c", gene_start, gene_end, f"{Gene}c", f"{Gene},{Partner}", "simMerged"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
                    out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in sim_merged_B.items():
            if len(exon_list)>0:
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3]
                for j, exon in enumerate(exon_list):
                    merged = exon.copy()
                    Partner = getMate(Gene, -1)
                    exon += [Gene, gene_start, gene_end, f"{Partner}c", f"{Partner},{Gene}", "simMerged"]
                    merged += [f"{Partner}c", gene_start, gene_end, f"{Partner}c", f"{Partner},{Gene}", "simMerged"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
                    out.write("\t".join([str(x) for x in merged]) + "\n")
        for Gene, exon_list in unchanged.items():
            if len(exon_list)>0:
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3]
                for j, exon in enumerate(exon_list):
                    exon += [Gene, gene_start, gene_end, Gene, f"{Gene},{Gene}", "unchanged"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
        for Gene, exon_list in real.items():
            if len(exon_list)>0:
                gene_start = exon_list[0][2]
                gene_end = exon_list[-1][3]
                for j, exon in enumerate(exon_list):
                    exon += [Gene, gene_start, gene_end, real_info[Gene][4], real_info[Gene][5], "real"]
                    out.write("\t".join([str(x) for x in exon]) + "\n")
    return()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "parse split-gene candidates into bed format for bedtools intersect with longest transcript gff.")
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='reciprocal split-gene candidates from the syntenic homology pipeline')
    parser.add_argument('-o', type=str, metavar='output_dir', required=True, help='full path')
    parser.add_argument('-a', type=str, metavar='annotations', required=True, help='comma-separated list (no spaces) of the two annotations (with genes and exons) that correspond to the two annotations compared in the input file')
    parser.add_argument('-n', type=str, metavar='sample_ids', required=True, help='comma-separated list (no spaces) of the genotype names that correspond to the two annotations (same order).  This should match naming convention used for naming samples in RNAseq output')
    parser.add_argument('-t', type=float, metavar='tandem_dup_threshold', default = 0.1, help='cds size - the number of base pairs both genes overlap over the cds / the size of the cds')
    parser.add_argument('-e', type=int, metavar='minimum_exons', default = 4, help="minimum number of exons required to make a simulated split gene")
    parser.add_argument('-S', type=float, metavar='split_proportion', required=False, default=0.2, help="proportion of genes for which we want to make simulated splits")
    parser.add_argument('-M', type=float, metavar='merged_proportion', required=False, default=0.3, help="proportion of genes for which we want to make simulated merge")
    # parser.add_argument('-s', type=float, metavar='syntenic_genes', required=True, help="File containing syntenic genes for B and P. One gene name per line")
    parser.add_argument('-v', action="store_true")

    args = parser.parse_args()

    if not os.path.exists(args.o): os.makedirs(args.o)


    #PRINT TOTAL NUMBER OF CANDIDATES FOR EACH RUN AT END

    #This needs to return a list holding all of the information for split genes that we intend to print later
    real_parsed = parseFile(args.i, args.t)

    annotations = args.a.split(",")
    names = args.n.split(",")    

    for i, annot in enumerate(annotations):
        outfile = f"{args.o}/{names[i]}splits_{''.join(names)}_td{args.t}_m{args.e}"
        sim_splits, sim_merged_A, sim_merged_B, unchanged, real = getExons(annot, real_parsed, args.S, args.M)
        writeSim(sim_splits, sim_merged_A, sim_merged_B, unchanged, real, real_parsed, outfile, args.e)

