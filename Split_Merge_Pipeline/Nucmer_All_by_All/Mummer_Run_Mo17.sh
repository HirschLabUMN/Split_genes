#!/bin/bash

#PBS -l walltime=24:00:00,nodes=1:ppn=2,pmem=50000mb
#PBS -V
#PBS -N MUMMER_RUN
#PBS -M mich0391@umn.edu
#PBS -m abe
#PBS -r n

#load GNU parallel

cd /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/Mummer_Results
#go to where the parallelization command is stored
bash Mo17_nucmer_commands.sh
