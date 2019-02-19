import subprocess
import os
import sys
import getopt


def usage():
    print("""\n
        This is the usage function:
        python3 Preprocess/Get_L_fasta.py -i <input_file> -f <fasta_file> -o <output_file>
            -i or --input_file   : Name of your table input file
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "i:f:o:h", ["input_file=",
                                           "fasta_file="
                                           "output_file="
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
    elif opt in ("-f", "--fasta_file"):
        fasta_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

data = open(infile, "r")


# Make calling a subprocess easier
def subprocess_cmd(command):
    subprocess.call(command,
                    stdout=subprocess.PIPE,
                    shell=True)


for input_line in data:
    input_line = input_line.strip()
    results = subprocess_cmd("sed -n '/" +
                             input_line +
                             "/,/^>/p' " +
                             fasta_file +
                             " | head -n -1 >> " +
                             output_file
                             )
