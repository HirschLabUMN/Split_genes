import sys
import getopt


def usage():
    print("""\n
        This is the usage function:
        python3 GFF_process.py  -i <input_gff>
            -i or --input_file   : Name of your table input file
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "i:h", ["input_file=",
                                       "help"]
                               )

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--input_file"):
        infile = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

# set dummy variables
Chr = 0
gene_count = 1
data = open(infile, "r")
for row in data:
    row = row.strip()
    if row.startswith("#"):
        continue
    else:
        fields = row.split("\t")
        # we only care about the genes in the gff
        if fields[2] == "gene":
            # reset counter if we find a gene on a new chr
            if Chr != fields[0]:
                Chr = fields[0]
                gene_count = 1
            fields.append(gene_count)
            print(*fields, sep="\t")
            gene_count += 1
