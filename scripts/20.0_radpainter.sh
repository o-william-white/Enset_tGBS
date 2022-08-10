#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_radpainter

module load gsl
module add R/3.6.1
export PATH=/data/home/mpx469/software/fineRADstructure/:$PATH

mkdir -p radpainter

# cp radpainter output
ln -sf ${PWD}/populations/populations_80_all_snps_clean/populations.haps.radpainter ${PWD}/radpainter/

# clone correct input
Rscript additional_scripts/clone_correct_radpainter_input.R

# radpainter
RADpainter paint radpainter/populations.haps.clone.correct.radpainter

# finestructure
finestructure -x 100000 -y 100000 -z 1000 radpainter/populations.haps.clone.correct_chunks.out radpainter/populations.haps.clone.correct_chunks.mcmc.xml
finestructure -m T -x 10000 radpainter/populations.haps.clone.correct_chunks.out radpainter/populations.haps.clone.correct_chunks.mcmc.xml radpainter/populations.haps.clone.correct_chunks.mcmcTree.xml

# plot
Rscript additional_scripts/plot_coancestry_radpainter.R

echo done

