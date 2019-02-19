input = open("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/"
             "Blast/Maize/CDS_results/BTX_623_LongestTS.fasta",
             "r")

for line in input:
    line = line.strip()
    if line.startswith(">"):
        fields = line.split(" ")
        fields[0] = fields[0].rsplit(".", 1)[0]
        fields = " ".join(fields)
        print(fields)
    else:
        print(line)
