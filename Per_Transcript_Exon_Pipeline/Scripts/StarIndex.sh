# Make your DIRs to house your STAR indexing files
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_B73
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_PH207
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_W22
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_BTX_623
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_B73_LT
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_PH207_LT
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_W22_LT
mkrir -p /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_BTX_623_LT

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_B73 \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa \
--sjdbGTFfile /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_PH207 \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa \
--sjdbGTFfile /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_W22 \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta \
--sjdbGTFfile /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/W22/zea_maysw22_core_fixchr.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_BTX_623 \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Sorghum_bicolor/btx623/Sbicolor_313_v3.0.fa \
--sjdbGTFfile /panfs/roc/groups/6/maize/shared/databases/genomes/Sorghum_bicolor/btx623/Sbicolor_313_v3.1.gene_exons.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

# Longest Transcript indexing
STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_B73_LT \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa \
--sjdbGTFfile /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Zea_mays.AGPv4.33_Longest_transcript_Exons.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_PH207_LT \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa \
--sjdbGTFfile /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/ZmaysPH207_443_v1.1.gene_exons_Longest_Transcript_Exon.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_W22_LT \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta \
--sjdbGTFfile /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/zea_maysw22_core_fixchr_Longest_Transcript_Exon.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50

STAR --runThreadN 24 \
--limitGenomeGenerateRAM=55000000000 \
--runMode genomeGenerate \
--genomeDir /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Index/STAR_50_BTX_623_LT \
--genomeFastaFiles /panfs/roc/groups/6/maize/shared/databases/genomes/Sorghum_bicolor/btx623/Sbicolor_313_v3.0.fa \
--sjdbGTFfile /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Sbicolor_313_v3.1.gene_exons_Longest_Transcript_Exon.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 50
