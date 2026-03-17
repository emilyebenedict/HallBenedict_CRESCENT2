#!/bin/bash

#===============================================================================
# File Name    : panaroo_bakta.sh
# Description  : Runs the panaroo pangenome tool on genomes annotated by bakta
# Usage        : sbatch panaroo_bakta.sh
#===============================================================================

#SBATCH --job-name=internal_panaroo
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --output=slurm_out/panaroo/x_panaroo_%A.out
#SBATCH --error=slurm_out/panaroo/y_panaroo_%A.err

eval $( spack load --sh py-panaroo@1.2.10 )

basedir="$PWD"
indir="${basedir}/d07_bakta/internal_gffs"
outdir="${basedir}/d16_panaroo_bakta/all_internal"

mkdir -p ${outdir}

set -x
time panaroo \
        -i ${indir}/*.gff3 \
        -o ${outdir} \
        --clean-mode strict \
        --core_threshold 0.95 \
        -a core \
        --aligner mafft \
        --remove-invalid-gene \
        -t ${SLURM_CPUS_PER_TASK}

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi
