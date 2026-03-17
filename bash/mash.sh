#!/bin/bash
#===============================================================================
#
# File Name    : s00_mash.sh
# Description  : This script will screen a genome with MASH
# Usage        : sbatch s00_mash.sh
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=mash
#SBATCH --array=1-3
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --output=slurm_out/mash/z_mashScreen_%a.out
#SBATCH --error=slurm_out/mash/z_mashScreen_%a.out


eval $( spack load --sh mash@2.3 )

basedir="$PWD"
refdir="/ref/gdlab/data/mash_sketches"
indir="${basedir}/d03_spades"
outdir="${basedir}/d05_mash_reseq"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/reseq_sr_list.txt`

set -x
mash screen -w -p ${SLURM_CPUS_PER_TASK} ${refdir}/refseq.genomes.k21s1000_210525.msh ${indir}/${sample}/scaffolds.fasta  > ${outdir}/${sample}_tmp.tab
RC=$?
set +x

sort -gr ${outdir}/${sample}_tmp.tab > ${outdir}/${sample}.tab
rm ${outdir}/${sample}_tmp.tab

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
