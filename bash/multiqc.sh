#!/bin/bash

#===============================================================================
# Name         : s02_multiqc.sh
# Description  : Consolidates output from s01_fastqc.sh
# Usage        : sbatch s02_multiqc.sh
#===============================================================================

#Submission script for HTCF
#SBATCH --job-name=multiqc
#SBATCH --mem=4G
#SBATCH --output=slurm_out/multiqc/z_multiqc_%A.out
#SBATCH --error=slurm_out/multiqc/z_multiqc_%A.out

eval $( spack load --sh py-multiqc@1.27 )

basedir="$PWD"
indir="${basedir}/d02_read_qc"
rawin="${indir}"
cleanin="${indir}/clean_reads"
rawout="${indir}/multiqc"
cleanout="${indir}/multiqc/clean_reads"
mkdir -p ${rawout}
mkdir -p ${cleanout}

set -x

time multiqc ${rawin} -o ${rawout} -n raw_multiqc
time multiqc ${cleanin} -o ${cleanout} -n clean_multiqc

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi
