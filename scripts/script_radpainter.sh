#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_radpainter

module load gsl
module add R
export PATH=/data/home/mpx469/software/fineRADstructure/:$PATH

# cp radpainter output
cp /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_all_snps_blacklist_output/populations.haps.radpainter .

# clone correct input
Rscript clone_correct_radpainter_input.R

# radpainter
RADpainter paint populations.haps.clone.correct.radpainter

# finestructure
finestructure -x 100000 -y 100000 -z 1000 populations.haps.clone.correct_chunks.out populations.haps.clone.correct_chunks.mcmc.xml
finestructure -m T -x 10000 populations.haps.clone.correct_chunks.out populations.haps.clone.correct_chunks.mcmc.xml populations.haps.clone.correct_chunks.mcmcTree.xml

echo done

