#!/bin/bash

#===============================================================================
# Name         : s06_checkm.sh
# Description  : This script will run checkm on a spades output dir
# Usage        : sbatch s06_checkm.sh
#===============================================================================

#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=checkm
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --output=slurm_out/checkm/z_checkm_%A.out
#SBATCH --error=slurm_out/checkm/z_checkm_%A.out

eval $( spack load --sh py-checkm-genome )

basedir="$PWD"
indir="${basedir}/d03_spades"
tmpdir="${basedir}/tmp_checkM_reseq"
outdir="${basedir}/d04_checkm_reseq"

mkdir -p ${outdir}

# copy scaffold files from unicycler to a temporary directory
# (checkM requires all your files being in a single folder)
# IMPORTANT: the find command is meant to find the output dirs from unicycler.
#mkdir -p ${tmpdir}
#for dir in `find ${indir} -type d -maxdepth 1 -mindepth 1`; do
#  cp $dir/scaffolds.fasta $tmpdir/$(basename $dir).fasta
#done

set -x
time checkm lineage_wf -f ${outdir}/CRESCENT_reseq_checkM_output.txt \
    -t ${SLURM_CPUS_PER_TASK} \
    -x fasta \
    --tab_table \
    ${tmpdir} ${outdir}
RC=$?
set +x


if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi
