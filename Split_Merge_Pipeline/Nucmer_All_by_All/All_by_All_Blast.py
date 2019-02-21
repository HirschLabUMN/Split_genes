import os
import sys
import tempfile
import subprocess
import pandas as pd
import getopt


# Listing all of the functions
def usage():
    print("""\n
        This is the usage function:
        All_by_All_Blast.py
            -q or --query_gff          : Path to qeury gff file (genes only)
            -s or --subject_gff        : Path to subject gff file (genes only)
            -b or --blast_DB           : Path to subject blast DB of longest
                                         CDS transcript
            -l or --longest_CDS        : path to fasta file of subjects longest
                                         representative transcript sequence
            -o or --output_file        : output file name
            -n or --nucmer_coords      : Path to nucmer .coords file between
                                         subject and query
            -g or --subject_char_split : Text delimiter of split gene character

            -h or --help               : help command
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "q:s:b:l:o:n:g:h", ["query_gff=",
                                                   "subject_gff=",
                                                   "blast_DB=",
                                                   "longest_CDS=",
                                                   "output_file="
                                                   "nucmer_coords="
                                                   "subject_char_split=",
                                                   ]
                               )

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-q", "--query_gff"):
        query_gff = arg
    elif opt in ("-s", "--subject_gff"):
        Subject_gff_file = arg
    elif opt in ("-b", "--blast_DB"):
        DB = arg
    elif opt in ("-l", "--longest_CDS"):
        CDS = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-n", "--nucmer_coords"):
        nucmer_coords = arg
    elif opt in ("-g", "--subject_char_split"):
        gene_char_split = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


def subprocess_cmd(command):
    subprocess.call(command,
                    stdout=subprocess.PIPE,
                    shell=True)


def cds_blast(gene, cds_fasta, DB):
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
                   cds_fasta +
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

    # Check if query file had anything
    if os.stat((base_split + 'query.fa')).st_size == 0:
        blast_genes = "CDS sequence not found"
        # remove old files
        subprocess_cmd('rm ' + base_split + 'query.fa')
        subprocess_cmd('rm ' + gene_file.name)
        gene_file.close()
        return(blast_genes)
    else:
        # remove old files
        subprocess_cmd('rm ' + base_split + 'query.fa')
        subprocess_cmd('rm ' + gene_file.name)
        gene_file.close()
        return (blast_genes)


def gene_grep(gene, gff):
    process = subprocess.check_output('grep -w ' +
                                      gene +
                                      ' ' +
                                      gff,
                                      shell=True
                                      ).decode('utf-8')
    process = process.strip()
    fields = process.split("\t")
    chrom = fields[0]
    start = fields[3]
    stop = fields[4]
    # pull out the gene name
    gene = fields[8].split(";")[0].split("=")[1]
    return(gene, chrom, start, stop)


def find_nearest_gene(blast_genes_file, syntenty_metadata):
    print("hi")


# End of functions, start script
# Get gene names from Gff file
Query_gff_file = open(query_gff, "r")

# Open results file
OutFile = open(output_file, "w")
# open the coordinates dataframe
df = pd.read_csv(nucmer_coords,
                 sep="\t",
                 index_col=None
                 )

# chromosome naming can confuse python since it can start as int then go to str
df.TAGS = df.TAGS.astype(str)
df.TAGS2 = df.TAGS2.astype(str)

File_Header = print("Query_gene",
                    "Q_Chr",
                    "Q_start",
                    "Q_stop",
                    "Q_Syn_Chr",
                    "Q_Syn_start",
                    "Q_Syn_stop",
                    "S_Syn_Chr",
                    "S_Syn_start",
                    "S_Syn_stop",
                    "#_genes_in_syntentic",
                    "Sytentic_genes",
                    "#_adjacent_genes",
                    "adjacent_genes",
                    "gene_type",
                    "CD2S_low",
                    "CD2S_high",
                    sep="\t",
                    file=OutFile
                    )

Gff_Header = Query_gff_file.readline()
for QueryGene in Query_gff_file:
    # QueryGene = Query_gff_file.readline()
    Q_fields = QueryGene.split()
    Q_chr = Q_fields[0]
    Q_start = Q_fields[3]
    Q_stop = Q_fields[4]
    Q_gene = Q_fields[8].split(";")[0].split("=")[1]
    # This is just for test cases
    # Q_chr = '4'
    # Q_start = 50276
    # Q_stop = 52356
    # Q_gene = 'Zm00001d027233'
    Q_gene_meta = [Q_gene, Q_chr, Q_start, Q_stop]

    # Call the function from above
    Blast_genes = set(cds_blast(Q_gene, CDS, DB))

    # Set vars as NA unless syntenic region not found
    CD2S_low = "NA"
    CD2S_high = "NA"

    # check to see if query gene was not found
    if "C" in Blast_genes:
        print(*Q_gene_meta,
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "CDS_sequence_for_query_gene_not_found",
              "NA",
              "NA",
              sep="\t",
              file=OutFile
              )
        Blast_genes = []
        continue
    else:
        # Get Chr, start and stop info for genes
        Blast_metadata = []
        for B_gene in Blast_genes:
            gene = B_gene.split("_")[0]
            B_gene_info = gene_grep(gene, Subject_gff_file)
            Blast_metadata.append(list(B_gene_info))

        # find region of syntenty for the query gene
        try:
            Dframe_index = df.index[(df["TAGS"] == Q_chr) &
                                    (df['S1'] <= int(Q_stop))].tolist()[-1]
            Syn = df.iloc[(Dframe_index)].tolist()
            Two_sided = True

        # This pipeline will fail if the first gene in a Chr is not in a
            # Syntenic region, this is why we set this index exception
        # If this fails then we can only do a one sided one
        except IndexError as e:
            print("type error: ",
                  str(e),
                  "for gene",
                  Q_gene,
                  "at begining of chr ",
                  "not in syntenic region")
            Two_sided = False
            continue

        Q_Syn_region_start = int(Syn[0])
        Q_Syn_region_stop = int(Syn[1])
        Q_Syn_region_chrom = Syn[11]
        Syn_region_start = int(Syn[2])
        Syn_region_stop = int(Syn[3])
        Syn_region_chrom = Syn[12]

        Syn_Meta = [Q_Syn_region_chrom,
                    Q_Syn_region_start,
                    Q_Syn_region_stop,
                    Syn_region_chrom,
                    Syn_region_start,
                    Syn_region_stop,
                    ]

        Syn_genes = []
        for result in Blast_metadata:
            gene = result[0]
            chrom = result[1]
            start = int(result[2])
            stop = int(result[3])
            if chrom == Syn_region_chrom:
                if (Syn_region_start) <= start <= (Syn_region_stop) or \
                   (Syn_region_start) <= stop <= (Syn_region_stop):
                    Syn_genes.append(result)
        #    print(len(Syn_genes))
        # print(Syn_genes)
        Syn_genes.sort(key=lambda x: x[0])

        if len(Syn_genes) > 1:
            adjacent_genes = []
            coordinates = []
            # print(Syn_genes)
            for idx, gene in enumerate(Syn_genes):
                coordinates.append(gene[0].split(gene_char_split)[1])
                if idx > 0:
                    if int(coordinates[idx]) - int(coordinates[idx-1]) < 3:
                        if Syn_genes[idx-1] not in adjacent_genes:
                            adjacent_genes.append(Syn_genes[idx-1])
                        adjacent_genes.append(Syn_genes[idx])
            if len(adjacent_genes) == 0:
                gene_type = "non-adjacent_syntenic"
                adjacent_genes_len = 0
            else:
                gene_type = "adjacent_genes_syntenic"
                adjacent_genes_len = len(adjacent_genes)

        elif len(Syn_genes) == 1:
            gene_type = "one_to_one_mapping"
            adjacent_genes = "NA"
            adjacent_genes_len = 0

        # Check to see if blast genes are near syntenic region
        else:
            # Set arbitrarily large values as a place holder
            CD2S_low = 1000000000000000
            CD2S_high = 1000000000000000
            # Check to see if you need to look on both sides of syntenic region
            if Two_sided:
                for index, B_gene in enumerate(Blast_metadata):
                    Before_candidate_list = []
                    After_candidate_list = []
                    if B_gene[1] == Syn_region_chrom:
                        before = Syn_region_start - int(B_gene[3])
                        after = int(B_gene[2]) - Syn_region_stop
                        if (before > 0) & (before < CD2S_low):
                            CD2S_low = before
                        elif (after > 0) & (after < CD2S_high):
                            CD2S_high = after
                        else:
                            continue
            else:
                for index, B_gene in enumerate(Blast_metadata):
                    After_candidate_list = []
                    if B_gene[1] == Syn_region_chrom:
                        after = int(B_gene[2]) - Syn_region_stop
                        if (after > 0) & (after < CD2S_high):
                            CD2S_high = after
                        else:
                            continue
            if CD2S_low == 1000000000000000:
                CD2S_low = "NA no CD2S low"
            if CD2S_high == 1000000000000000:
                CD2S_high = "NA no CD2S high"

            gene_type = "non-matching_syntenic"
            adjacent_genes = "NA"
            adjacent_genes_len = 0

        # format the Syn gene lists
        P_Syn_genes = []
        for gene in Syn_genes:
            P_Syn_genes.append(",".join(map(str, gene)))

        P_adjacent_genes = []
        for gene in adjacent_genes:
            P_adjacent_genes.append(",".join(map(str, gene)))

        print(*Q_gene_meta,
              *Syn_Meta,
              len(Syn_genes),
              ";".join(map(str, P_Syn_genes)),
              adjacent_genes_len,
              ";".join(map(str, P_adjacent_genes)),
              gene_type,
              CD2S_low,
              CD2S_high,
              sep="\t",
              file=OutFile
              )
        OutFile.flush()
# print(Q_gene, Q_start, Q_stop, Syn_genes)
OutFile.close()
