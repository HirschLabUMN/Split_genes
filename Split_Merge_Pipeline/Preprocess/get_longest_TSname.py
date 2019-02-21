input = open("/panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/"
             "Scripts/BTX623_Longest_transcript_Exon.gff",
             "r")

ts_name = ""
for line in input:
    if line.startswith("#"):
        continue
    else:
        line = line.strip()
        fields = line.split("\t")
        # only check genes
        if fields[2] == "gene":
            continue
        else:
            # This line will differ across GFF files
            new_ts = fields[8].split(";")[1].split("=")[1]
            # if our transript matches previous line, continue
            if new_ts == ts_name:
                continue
            else:
                # If not, print it and set as new val for ts_name
                print(new_ts)
                ts_name = new_ts
