#Manually declage the GFF Star indexes as well as the prefixes that we want to use
declare -a FASTA=("/home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Sorghum_bicolor/btx623/Sbicolor_313_v3.0.fa"
                  "/home/maize/shared/databases/genomes/Zea_mays/mo17_10chr_gbrowse_1August2018.fa"
                  "/home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/PHB47/Zea_mays_var_PHB47.mainGenome.fasta"
                  "/home/maize/shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta")

declare -a QUERY=("B73"
                  "BTX623"
                  "Mo17"
                  "PH207"
                  "PHB47"
                  "W22")

for j in {0..5}; do
    for i in {0..5}; do
      if [ "$i" -ne "$j" ]; then
        (echo nucmer --mum -c 1000 -p ${QUERY[$j]}_${QUERY[$i]}_c1000 ${FASTA[$j]} ${FASTA[$i]})
      fi
    done
done
