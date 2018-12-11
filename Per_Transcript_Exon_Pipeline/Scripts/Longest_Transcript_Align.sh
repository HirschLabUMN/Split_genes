
#!/bin/bash

#PBS -l walltime=24:00:00,nodes=15:ppn=24,pmem=2500mb
#PBS -V
#PBS -A hirschc1
#PBS -N jobLongTransAlign
#PBS -M mich0391@umn.edu
#PBS -m abe
#PBS -r n

#load GNU parallel
module load fastqc
module load samtools
module load parallel

#go to where the parallelization command is stored
cd /panfs/roc/scratch/michnoj0/SplitGenes
export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"
sort -u $PBS_NODEFILE > unique-nodelist.txt
parallel --jobs 1 --sshloginfile unique-nodelist.txt --workdir $PWD < /panfs/roc/groups/14/hirschc1/mich0391/Projects/SplitGenes/Scripts/ReadMapCommands_LT.txt
