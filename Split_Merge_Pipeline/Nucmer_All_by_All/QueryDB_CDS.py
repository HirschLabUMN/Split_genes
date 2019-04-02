import subprocess
import tempfile
import os
import sys
import getopt


def usage():
    print("""\n
        This is the usage function:
        python3 QueryDB.py
            -i or --input_file   : Name of your table input file
            -o or --outputfile   : Name of your results or output file
            -m or --merged_fgg   : path to merged CDS file
            -s or --split_gff    : path to split CDS file
            -d or --data_base    : path to balst data base location
            -h or --help         : help command
        \n""")

try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "i:o:m:s:d:h", ["input_file=",
                                               "output_file=",
                                               "merged_gff=",
                                               "split_gff=",
                                               "data_base="
                                               "help",
                                               ]
                               )

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--input_file"):
        infile = arg
    elif opt in ("-o", "--output_file"):
        outfile = arg
    elif opt in ("-m", "--merged_gff"):
        merge_gff = arg
    elif opt in ("-s", "--split_gff"):
        split_gff = arg
    elif opt in ("-d", "--data_base"):
        DB = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


# Set the hardcoded paths
data = open(infile, "r")
# results_file = open(outfile, "w")
# merge_gff = "B73_genes.gff3"
# merge_fasta = "/panfs/roc/groups/6/maize/shared/databases/genomes/" \
#               "Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa"
# split_gff = "W22_genes.gff3"
# split_fasta = "/home/maize/shared/databases/genomes/" \
#               "Zea_mays/W22/W22__Ver12.genome.normalized.fasta"


# merge_gene = "Zm00001d002216"
# splits = ["Zm00004b006245", "Zm00004b006246", "Zm00004b006244"]

# Write a header
print("Query_merged_gene",
      "Query_split_gene",
      "Num_split_genes",
      "State",
      "Call",
      "Output_call"
      "Num_genes_intersect",
      "Adjacent_genes_intersect",
      "Adjacent_genes",
      "Intersect_genes",
      "TandemDupCall",
      "TandemDupPer",
      "Adj_gene_meta",
      sep="\t")
sys.stdout.flush()
#results_file.flush()


# Make calling a subprocess easier
def subprocess_cmd(command):
    subprocess.call(command,
                    stdout=subprocess.PIPE,
                    shell=True)


def ortho_blast(gene, gff):
    # make a temp file
    gene_file = tempfile.NamedTemporaryFile(suffix="_temp.txt",
                                            prefix="gene_",
                                            mode='r+',
                                            delete=False,
                                            )
    # get the basename
    base = os.path.basename(gene_file.name)
    base_split = os.path.splitext(base)[0]
    # we just needed the unique name so we can close it
    # gene_file.close()

    #  sed command to pull gene from gff
    subprocess_cmd("sed -n '/" +
                   gene +
                   "/,/^>/p' " +
                   gff +
                   " | head -n -1 >> " +
                   base_split +
                   'query.fa')

    # blast to arabidopsis gene db
    subprocess_cmd('blastn ' +
                   '-db ' +
                   DB +
                   ' ' +
                   '-query ' +
                   base_split +
                   'query.fa ' +
                   '-out ' +
                   gene_file.name +
                   ' ' +
                   '-task blastn ' +
                   '-evalue 1e-4 ' +
                   '-outfmt 6')

    cds_fasta_file = open(base_split+"query.fa", "r")
    cds_length = 0
    for line in cds_fasta_file:
        line = line.strip()
        if line.startswith(">"):
            continue
        else:
            cds_length += len(line)
    cds_fasta_file.close()

    # open the blast results and return a set of genes
    blast_results_file = open(gene_file.name, "r")
    blast_genes = []
    pid = []
    qstart = []
    qend = []
    e_val = []
    for line in blast_results_file:
        line = line.strip()
        fields = line.split()
        blast_genes.append(fields[1])
        pid.append(fields[2])
        qstart.append(fields[6])
        qend.append(fields[7])
        e_val.append(fields[10])
    # make it a unique list
    # unique_blast_genes = set(blast_genes)
    # remove old files
    subprocess_cmd('rm ' + base_split + 'query.fa')
    subprocess_cmd('rm ' + gene_file.name)
    gene_file.close()

    return (blast_genes, pid, qstart, qend, e_val, cds_length)

for input_line in data:
    input_line = input_line.strip()
    input_fields = input_line.split("\t")
    merge_gene = input_fields[1]
    # splits = input_fields[2].split(",")
    splits = [i.split(",")[0] for i in input_fields[14].split(";")]
    call = input_fields[4]
    state = input_fields[5]

    M_Blast = ortho_blast(merge_gene,
                          merge_gff,
                          )
    merge_list = M_Blast[0]
    Blast_pid = M_Blast[1]
    Blast_Qstart = M_Blast[2]
    Blast_Qstop = M_Blast[3]
    Blast_e_val = M_Blast[4]
    cds_length = int(M_Blast[5])
    cds_range = range(1, (cds_length + 1))

    # make an overlap set
    split_list = []
    all_list = []
    for split in splits:
        split_results = ortho_blast(split,
                                    split_gff,
                                    )
        # take the intersect of sets and keep updating it
        split_list.append(split_results[0])
        all_list.append(split_results)

    all_genes = []
    all_pid = []
    all_Qstart = []
    all_Qstop = []
    all_e_val = []
    overlap_metadata = "NA"
    go_print = []
    t_dup_percent = "NA"
    for il in all_list:
        all_genes += il[0]
        all_pid += il[1]
        all_Qstart += il[2]
        all_Qstop += il[3]
        all_e_val += il[4]

# make sure that the split hit the same genes
    split_set = set.intersection(*[set(list) for list in split_list])
    merge_set = set(merge_list)

    # intersection of all three
    gene_overlap = list(merge_set.intersection(split_set))
    if len(gene_overlap) > 1:
        Blast_metadata = []
        for idx, B_gene in enumerate(all_genes):
            idxs = [i for i, x in enumerate(all_genes) if x == B_gene]
            if len(idxs) > 1:
                largest_val_index = 0
                largest_blast_size = 0
                smallest_e_val = float(10000)
                # Filter for multiple blasts to same gene by eval then len
                for idx in idxs:
                    blast_size = (int(all_Qstop[idx]) -
                                  int(all_Qstart[idx]))
                    current_e_val = float(all_e_val[idx])
                    if current_e_val < smallest_e_val:
                        largest_val_index = idx
                        smallest_e_val = current_e_val
                    # in the off chance of equal e-vals
                    elif current_e_val == smallest_e_val:
                        if blast_size > largest_blast_size:
                            largest_blast_size = blast_size
                            largest_val_index = idx

                B_gene_info = (all_genes[largest_val_index],
                               all_pid[largest_val_index],
                               all_Qstart[largest_val_index],
                               all_Qstop[largest_val_index],
                               all_e_val[largest_val_index],
                               )
            Blast_metadata.append(list(B_gene_info))

        # OK, this is where things get tricky, we need to go back to the blast
        # metadata, pull out all the information relating to genes,
        # Then filter for multiple blast hits. lets get this show on the road
        # This method is going to be super slow I don't have the time right now
        # Since we are calling this split, we expect two distinct genes
        # so we can append everything to one list and filter for duplicate
        # gene blasts later
        # I apologize to whoever reads this code
        adjacent_genes = []
        atg_chrs = []
        coordinates = []
        gene_overlap.sort()

        for idx, gene in enumerate(gene_overlap):
            atg_chrs.append(gene.split("G")[0])
            semicoord = gene.split("G")[1]
            coordinates.append(semicoord.split(".")[0])
            if idx > 0 and atg_chrs[idx] == atg_chrs[idx-1]:
                if int(coordinates[idx]) - int(coordinates[idx-1]) < 50:
                    if gene_overlap[idx-1] not in adjacent_genes:
                        adjacent_genes.append(gene_overlap[idx-1])
                    adjacent_genes.append(gene_overlap[idx])

            if len(adjacent_genes) == 0:
                padjacent_genes = "NA"
                # output_call = "NOT CALLABLE"

        go = []
        lst2 = [item[0] for item in Blast_metadata]
        for p in adjacent_genes:
            index = lst2.index(p)
            go.append(Blast_metadata[index])

        go_len = len(go)
        ranges = []
        s_cds_range = set(cds_range)
        if go_len > 1:

            for gene in go:
                go_print.append(",".join(gene))
                blast_range = list(range(int(gene[2]),
                                   (1+int(gene[3]))))
                ranges.append(blast_range)

        if go_len == 2:
            # use sets because we are useing buil in function
            R1 = set(ranges[0])
            R2 = set(ranges[1])
            two_overlap = len(set.intersection(s_cds_range,
                                               R1,
                                               R2
                                               ))
            one_overlap = (len(set.intersection(s_cds_range, R1)) +
                           len(set.intersection(s_cds_range, R2)) -
                           (two_overlap*2
                            ))
            zero_overlap = cds_length - one_overlap - two_overlap
            three_overlap = "NA"

            # tandem dup percentage calc
            t_dup_percent = (two_overlap) / (two_overlap +
                                             one_overlap)
            overlap_list = [str(cds_length),
                            str(zero_overlap),
                            str(one_overlap),
                            str(two_overlap),
                            str(three_overlap)]
            overlap_metadata = ":".join(overlap_list)

        elif go_len == 3:
            R1 = set(ranges[0])
            R2 = set(ranges[1])
            R3 = set(ranges[2])
            three_overlap = len(set.intersection(s_cds_range,
                                                 R1,
                                                 R2,
                                                 R3,))
            two_overlap = (len(set.intersection(s_cds_range,
                                                R1,
                                                R2,
                                                )) +
                           len(set.intersection(s_cds_range,
                                                R2,
                                                R3)) +
                           len(set.intersection(s_cds_range,
                                                R1,
                                                R3)) -
                           (three_overlap*3))

            one_overlap = (len(set.intersection(s_cds_range, R1)) +
                           len(set.intersection(s_cds_range, R2)) +
                           len(set.intersection(s_cds_range, R3)) -
                           (three_overlap*3) - (two_overlap * 2))
            zero_overlap = (cds_length -
                            one_overlap -
                            two_overlap -
                            three_overlap)

            # tandem dup percentage calc
            t_dup_percent = (three_overlap) / (three_overlap +
                                               two_overlap +
                                               one_overlap)

            overlap_list = [str(cds_length),
                            str(zero_overlap),
                            str(one_overlap),
                            str(two_overlap),
                            str(three_overlap)]
            overlap_metadata = ":".join(overlap_list)

        elif go_len > 3:
            overlap_list = [str(cds_length),
                            "too many combinations",
                            "NA",
                            "NA"]
            overlap_metadata = ":".join(overlap_list)
            t_dup_percent = "too many combinations"

        padjacent_genes = ",".join(adjacent_genes)
        ladjacent_genes = len(adjacent_genes)

        output_call = "SPLIT"
        if ladjacent_genes == 0:
            padjacent_genes = "NA"
            output_call = "NOT CALLABLE"

    elif len(gene_overlap) == 1:
        ladjacent_genes = 0
        padjacent_genes = "NA"
        output_call = "MERGED"

    else:
        ladjacent_genes = 0
        padjacent_genes = "NA"
        output_call = "NOT CALLABLE"

    pgene_overlap = ",".join(gene_overlap)
    pgene_splits = ",".join(splits)

    print(merge_gene,
          pgene_splits,
          len(splits),
          state,
          call,
          output_call,
          len(gene_overlap),
          ladjacent_genes,
          padjacent_genes,
          pgene_overlap,
          overlap_metadata,
          t_dup_percent,
          (";".join(go_print)),
          sep="\t")
    sys.stdout.flush()
#    results_file.flush()


# results_file.close()
data.close()
