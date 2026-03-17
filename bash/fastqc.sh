#!/bin/bash
#===============================================================================
# Name         : s01_fastqc.sh
# Description  : Runs basic QC on raw and trimmed illumina reads
# Usage        : sbatch s00_fastQC.sh
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=fastqc
#SBATCH --array=1-3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm_out/fastqc/z_fastqc_%a_%A.out
#SBATCH --error=slurm_out/fastqc/z_fastqc_%a_%A.out

# load latest version of fastqc on spack (fastqc@0.11.9)
eval $( spack load --sh /3ng5m3l )

basedir="$PWD"
rawin="/lts/gdlab/users/current/ebenedict/crescent/d00_reseq_raw"
cleanin="${basedir}/d01_clean_reads"
rawout="${basedir}/d02_read_qc"
cleanout="${basedir}/d02_read_qc/clean_reads"
mkdir -p ${rawout}
mkdir -p ${cleanout}
export JAVA_ARGS="-Xmx8000M"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/reseq_sr_list.txt`

set -x

time fastqc ${rawin}/${sample}_R1_001.fastq.gz \
            ${rawin}/${sample}_R2_001.fastq.gz \
            -o ${rawout} \
            -t ${SLURM_CPUS_PER_TASK}

time fastqc ${cleanin}/${sample}_FW_clean.fastq.gz \
            ${cleanin}/${sample}_RV_clean.fastq.gz \
            ${cleanin}/${sample}_UP_clean.fastq.gz \
            -o ${cleanout} \
            -t ${SLURM_CPUS_PER_TASK}

RC=$?
set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
