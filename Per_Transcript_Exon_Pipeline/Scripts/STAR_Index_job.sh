#!/bin/bash

#PBS -l walltime=1:30:00,nodes=1:ppn=24,pmem=2500mb
#PBS -V
#PBS -A hirschc1
#PBS -N STAR_index_LT
#PBS -M mich0391@umn.edu
#PBS -m abe
#PBS -r n

cd /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes

bash /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/StarIndex.sh
