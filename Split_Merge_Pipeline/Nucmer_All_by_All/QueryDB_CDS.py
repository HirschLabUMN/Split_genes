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
            -d or --data_base    : path to blast data base location
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
results_file = open(outfile, "w")

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
      sep="\t",
      file=results_file)
results_file.flush()


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

    # open the blast results and return a set of genes
    blast_results_file = open(gene_file.name, "r")
    blast_genes = []
    for line in blast_results_file:
        line = line.strip()
        fields = line.split()
        blast_genes.append(fields[1])
    # make it a unique list
    # unique_blast_genes = set(blast_genes)
    # remove old files
    subprocess_cmd('rm ' + base_split + 'query.fa')
    subprocess_cmd('rm ' + gene_file.name)
    gene_file.close()

    return (blast_genes)


for input_line in data:
    input_line = input_line.strip()
    input_fields = input_line.split("\t")
    merge_gene = input_fields[0]
    # splits = input_fields[2].split(",")
    splits = [i.split(",")[0] for i in input_fields[13].split(";")]
    call = input_fields[3]
    state = input_fields[4]

    merge_list = ortho_blast(merge_gene,
                             merge_gff,
                             )

    # make an overlap set
    split_list = []
    for split in splits:
        split_results = ortho_blast(split,
                                    split_gff,
                                    )
        # take the intersect of sets and keep updating it
        split_list.append(split_results)

    split_set = set.intersection(*[set(list) for list in split_list])
    merge_set = set(merge_list)

    gene_overlap = list(merge_set.intersection(split_set))
    if len(gene_overlap) > 1:
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
                output_call = "NOT CALLABLE"

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
          sep="\t",
          file=results_file)
    results_file.flush()


results_file.close()
data.close()
