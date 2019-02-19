
#!/bin/bash

#PBS -l walltime=36:00:00,nodes=1:ppn=4,pmem=2500mb
#PBS -V
#PBS -N adopDB_Q
#PBS -M mich0391@umn.edu
#PBS -m abe
#PBS -r n

#load GNU parallel
module load bedtools
module load bcftools
module load parallel

#go to where the parallelization command is stored
cd /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Maize
parallel --jobs 6  < pycommands.txt
