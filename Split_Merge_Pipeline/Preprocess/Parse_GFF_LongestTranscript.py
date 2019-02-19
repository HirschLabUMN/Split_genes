"""This sctip will take a gff file as input and parce out the longest
   transcripts  in gff format for formatting for exon mapping"""

import sys
import getopt
import operator


def usage():
    print("""\n
        This is the usage function:
        python3 Parse_GFF_LongestTranscript.py -c <COB> -g <GOnt>
            -g or gff_file The name of your co-expression network
            -o or output_file  The name of your GO network object
            -i or id_field the value of field in gff representative of
                           transcript id (typically "Parent" or "ID")
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:i:h", ["gff_file=",
                                                         "output_file=",
                                                         "id_field="
                                                         "help"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-g", "--gff_file"):
        gff_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-i", "--id_field"):
        Transcript_ID = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

# open file for reading, writing and ID len
in_file = open(gff_file, "r")
out_file = open(output_file, "w")
T_ID_Length = len(Transcript_ID)

# initialize the first round  of dicts even though it will get looped over
Transcript_Compare = dict()
Gff_line_exon = dict()
for row in in_file:
    line = row.rstrip()
    # Keep all of the header files
    if line.startswith("#"):
        print(row, file=out_file, end="")
    else:
        # We are just looking for information then printing the entire row
        # so we are not going to strip newline
        fields = line.split("\t")
        # If we find a gene, this signifies the step to find the longest
        # transcript for the prevoius gene
        if fields[2] == "gene":
            # This is to prevent erros on first loop with
            if bool(Transcript_Compare):
                # Return key with larges value
                Longest_transcript = max(Transcript_Compare,
                                         key=Transcript_Compare.get)
                # Use key to pull out raw row informaiton for transcripts
                Exon_Lines = Gff_line_exon[Longest_transcript]
                print("".join(Exon_Lines), file=out_file, end="")

##########################################################################
                # # This is just for testing, uncomment if you want to see
                # # the length outputs for each transcript and the longest
                # for k, v in Transcript_Compare.items():
                #     print(k, v)
                # print("Longest=", Longest_transcript)
##########################################################################

            # print out the gene row information
            print(row, file=out_file, end="")
            # reinitialize/empty dicts for next gene
            Transcript_Compare = dict()
            Gff_line_exon = dict()

        elif fields[2] == "exon":
            # Get the size of the exon
            Size = int(fields[4]) - int(fields[3])
            # Pull out and split information from ID category
            ID_Fields = fields[8].split(";")
            for ID_Field_val in ID_Fields:
                # Look to see what item has your Transcript_ID keyword
                if Transcript_ID in ID_Field_val:
                    # pull out that name removing the prefix
                    Transcript = ID_Field_val[(T_ID_Length+1):]

            # Lood to see if that transcript is already in the dict
            if Transcript in Transcript_Compare.keys():
                Transcript_Compare[Transcript] += (Size)
                # Store the raw row information using the same key
                Gff_line_exon[Transcript].append(row)
            else:
                # if not namke a new one
                Transcript_Compare[Transcript] = Size
                Gff_line_exon[Transcript] = [row]

# Run for the last gene
if bool(Transcript_Compare):
    Longest_transcript = max(Transcript_Compare,
                             key=Transcript_Compare.get)
    # Use key to pull out raw row informaiton for transcripts
    Exon_Lines = Gff_line_exon[Longest_transcript]
    print("".join(Exon_Lines), file=out_file, end="")

    # print out the gene row information
    print(row, file=out_file, end="")

# close files for safety
in_file.close()
out_file.close()
