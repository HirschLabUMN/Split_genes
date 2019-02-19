import subprocess
import tempfile
import os
import sys
import getopt


def usage():
    print("""\n
        This is the usage function:
        python3 QueryDB.py -i <input_file> -o <outputfile> -s <query_gff>
                           -s <subject_gff>
            -i or --input_file   : Name of your table input file
            -o or --outputfile   : Name of your results or output file
            -q or --query_gff    : path to merged gff file
            -s or --subject_gff  : path to subject gff file
            -h or --help         : help command
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "i:o:q:s:f:e:d:h", ["input_file=",
                                                   "output_file=",
                                                   "query_gff=",
                                                   "subject_gff=",
                                                   "merged_fasta=",
                                                   "split_fasta=",
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
    elif opt in ("-q", "--query_gff"):
        query_gff = arg
    elif opt in ("-s", "--subject_gff"):
        subject_gff = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


# Set the hardcoded paths
data = open(infile, "r")
results_file = open(outfile, "w")
# merge_gff = "B73_genes.gff3"
# merge_fasta = "/panfs/roc/groups/6/maize/shared/databases/genomes/" \
#               "Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa"
# split_gff = "W22_genes.gff3"
# split_fasta = "/home/maize/shared/databases/genomes/" \
#               "Zea_mays/W22/W22__Ver12.genome.normalized.fasta"


# merge_gene = "Zm00001d002216"
# splits = ["Zm00004b006245", "Zm00004b006246", "Zm00004b006244"]

# Write a header
# print("Query_merged_gene",
#       "Query_split_gene",
#       "Num_split_genes",
#       "State",
#       "Call",
#       "Output_call"
#       "Num_genes_intersect",
#       "Adjacent_genes_intersect",
#       "Adjacent_genes",
#       "Intersect_genes",
#       sep="\t",
#       file=results_file)
# results_file.flush()


# Make calling a subprocess easier
def subprocess_cmd(command):
    subprocess.call(command,
                    stdout=subprocess.PIPE,
                    shell=True)


def gene_grep(gene, gff):
    # make a temp file
    gene_file = tempfile.NamedTemporaryFile(suffix="_temp.txt",
                                            prefix="gene",
                                            mode='r+',
                                            delete=False,
                                            )

    #  grep command to pull gene from gff
    subprocess_cmd('grep ' +
                   gene +
                   ' ' +
                   gff +
                   ' > ' +
                   # merge var and string for unique file name
                   gene_file.name)

    # open the blast results and return a set of genes
    results_file = open(gene_file.name, "r")
    gene_fields = results_file.readline().strip().split("\t")
    # pull out Chr, start, stop, and gene order # on Chr/scaffold
    metadata = [gene_fields[0],
                gene_fields[3],
                gene_fields[4],
                gene_fields[-1]]

    gene_file.close()
    subprocess_cmd('rm ' + gene_file.name)

    return (metadata)


for input_line in data:
    input_line = input_line.strip()
    input_fields = input_line.split("\t")
    query_gene = input_fields[0].split("_")[0]
    query_gene = query_gene.split("m.g")[0]
    query_gene = query_gene.rsplit('.', 1)[0]
    subject_gene = input_fields[1].split("_")[0]
    subject_gene = subject_gene.split("m.g")[0]
    subject_gene = subject_gene.rsplit('.', 1)[0]
    pval = input_fields[10]

    query_list = gene_grep(query_gene, query_gff)
    subject_list = gene_grep(subject_gene, subject_gff)

    print(query_list[0],
          query_gene,
          query_list[3],
          query_list[3],
          subject_list[0],
          subject_gene,
          subject_list[3],
          subject_list[3],
          pval,
          sep="\t",
          file=results_file)
    results_file.flush()
results_file.close()
