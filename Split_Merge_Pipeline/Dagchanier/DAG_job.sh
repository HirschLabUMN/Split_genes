#!/bin/bash

#PBS -l walltime=48:00:00,nodes=1:ppn=20,pmem=2500mb
#PBS -V
#PBS -N dagjobs
#PBS -M mich0391@umn.edu
#PBS -m abe
#PBS -r n

#load GNU parallel
module load parallel

#go to where the 
cd /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Blast/Dagchainer
parallel --jobs 20 < DAG_Commands.txt