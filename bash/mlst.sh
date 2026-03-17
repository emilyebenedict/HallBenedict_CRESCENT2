#!/bin/bash
#===============================================================================
# File Name    : mlst.sh
# Description  : Scan many genomes against PubMLST typing schemes
# Usage        : sbatch mlst.sh
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-0:00:00 # days-hh:mm:ss
#SBATCH --job-name=mlst@2.22.1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --output=slurm_out/mlst/z_mlst_%A.out
#SBATCH --error=slurm_out/mlst/z_mlst_%A.out

eval $( spack load --sh mlst )

basedir="$PWD"
indir="${basedir}/tmp_checkM"
outdir="${basedir}/d10_mlst"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/seq_list.txt`

set -x
time mlst --csv \
  --nopath \
  --threads ${SLURM_CPUS_PER_TASK} \
  ${indir}/*.fasta >> ${outdir}/CRESCENT_mlst_results.csv
RC=$?
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi
